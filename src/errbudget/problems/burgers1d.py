"""L2: 1D viscous Burgers -- one rung harder than linear Crank-Nicolson.

    u_t + u u_x = nu u_xx        on [0, 2pi), periodic,  u(x,0) = sin x

Why this rung:
  * Nonlinear, so each time step needs a Newton solve wrapped around a GMRES
    solve -- two *nested* tolerances, which is exactly the configuration
    external users get wrong.
  * A steep front forms at t = 1 (the inviscid breaking time), so observables
    are genuinely sensitive to spatial resolution.
  * Ground truth is still exact, via the Cole-Hopf transformation:
        theta_t = nu theta_xx ,  u = -2 nu theta_x / theta ,
        theta(x,0) = exp(-(1-cos x) / (2 nu)) ,
    and the heat equation is solved to machine precision in Fourier space.
    So we can validate against truth, not against a fine-grid surrogate.

Four parameters: N, steps, newton_tol, gmres_tol.
Two observables from the same run: energy and dissipation rate at T.  They are
not minimised by the same configuration -- which matters, because users pick
parameters once and then reuse them for every quantity they measure.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class Burgers1D:
    def __init__(self, nu: float = 0.05, T: float = 1.0, max_newton: int = 25):
        self.nu, self.T, self.max_newton = nu, T, max_newton

    # ------------------------------------------------------- exact solution
    def exact(self, t: float, Nf: int = 1 << 13, ny: int = 1 << 12, L: float = 12.0):
        """Cole-Hopf solution, evaluated in log space.

        The direct route -- diffuse theta(x,0) = exp(-(1-cos x)/2nu) spectrally
        and set u = -2 nu theta_x/theta -- loses all significance once nu is
        small, because theta spans exp(-1/nu) orders of magnitude while the
        transform's round-off is set by its largest entry.  At nu = 0.03 theta
        underflows into the FFT's own noise and u is returned as garbage.  The
        equivalent integral form

            u(x,t) = INT ((x-y)/t) e^{-G/2nu} dy  /  INT e^{-G/2nu} dy ,
            G(y; x,t) = (1 - cos y) + (x-y)^2 / 2t ,

        is stable for every nu once the exponentials are taken relative to
        min_y G, since only ratios of the weights enter.  The Gaussian factor
        confines the integrand to |x-y| = O(sqrt(nu t)), so a fixed window of
        half-width L is ample.
        """
        nu = self.nu
        x = 2 * np.pi * np.arange(Nf) / Nf
        s = np.linspace(-L, L, ny)                       # s = x - y
        u = np.empty(Nf)
        for a in range(0, Nf, 512):                      # chunked: Nf x ny is large
            xc = x[a:a + 512]
            y = xc[:, None] - s[None, :]
            G = (1.0 - np.cos(y)) + s[None, :] ** 2 / (2.0 * t)
            w = -G / (2.0 * nu)
            w -= w.max(axis=1, keepdims=True)            # log-sum-exp shift
            e = np.exp(w)
            u[a:a + 512] = (np.sum((s[None, :] / t) * e, axis=1)
                            / np.sum(e, axis=1))
        k = np.fft.fftfreq(Nf, d=1.0 / Nf)
        ux = np.real(np.fft.ifft(1j * k * np.fft.fft(u)))
        return x, u, ux

    def exact_observables(self, Nf: int = 1 << 13) -> dict:
        x, u, ux = self.exact(self.T, Nf)
        dx = 2 * np.pi / Nf
        return {"energy": 0.5 * float(np.sum(u**2)) * dx,
                "dissipation": self.nu * float(np.sum(ux**2)) * dx}

    # ------------------------------------------------------------ operators
    @staticmethod
    def _ops(N: int):
        h = 2 * np.pi / N
        o = np.ones(N)
        Dx = sp.diags([o, -o], [1, -1], shape=(N, N), format="lil") / (2 * h)
        Dx[0, -1], Dx[-1, 0] = -1 / (2 * h), 1 / (2 * h)
        Lx = sp.diags([o, -2 * o, o], [-1, 0, 1], shape=(N, N), format="lil") / h**2
        Lx[0, -1], Lx[-1, 0] = 1 / h**2, 1 / h**2
        return Dx.tocsr(), Lx.tocsr()

    def _F(self, u, Dx, Lx):
        return -u * (Dx @ u) + self.nu * (Lx @ u)

    def _J(self, u, Dx, Lx, dt, N):
        # F'(u) v = -(u v_x + v u_x) + nu v_xx
        Fp = -sp.diags(u) @ Dx - sp.diags(Dx @ u) + self.nu * Lx
        return (sp.identity(N, format="csr") - 0.5 * dt * Fp).tocsr()

    # ------------------------------------------------------------------ run
    def run(self, theta: dict):
        N = int(theta["N"])
        steps = int(theta["steps"])
        ntol = float(theta["newton_tol"])
        gtol = float(theta["gmres_tol"])
        dt = self.T / steps
        h = 2 * np.pi / N

        x = 2 * np.pi * np.arange(N) / N
        u = np.sin(x)
        Dx, Lx = self._ops(N)

        matvecs = newton_its = 0
        newton_fail = gmres_fail = 0      # steps / solves that missed their tolerance
        t0 = time.perf_counter()
        for _ in range(steps):
            un = u.copy()
            Fn = self._F(un, Dx, Lx)
            v = un.copy()
            scale = np.linalg.norm(un) + 1e-30
            for _ in range(self.max_newton):
                G = v - un - 0.5 * dt * (self._F(v, Dx, Lx) + Fn)
                newton_its += 1
                if np.linalg.norm(G) <= ntol * scale:
                    break
                J = self._J(v, Dx, Lx, dt, N)
                cnt = [0]
                cb = lambda _r: cnt.__setitem__(0, cnt[0] + 1)
                try:
                    dv, info = spla.gmres(J, -G, rtol=gtol, atol=0.0,
                                          restart=min(50, N), maxiter=400,
                                          callback=cb, callback_type="pr_norm")
                except TypeError:
                    dv, info = spla.gmres(J, -G, tol=gtol, atol=0.0,
                                          restart=min(50, N), maxiter=400,
                                          callback=cb, callback_type="pr_norm")
                gmres_fail += info > 0
                matvecs += max(cnt[0], 1)
                v = v + dv
            else:
                G = v - un - 0.5 * dt * (self._F(v, Dx, Lx) + Fn)
                newton_fail += np.linalg.norm(G) > ntol * scale
            u = v
        wall = time.perf_counter() - t0

        ux = Dx @ u
        obs = {"energy": 0.5 * float(np.sum(u**2)) * h,
               "dissipation": self.nu * float(np.sum(ux**2)) * h}
        # flop proxy: sparse matvecs dominate, each O(N)
        cost = float(N * (matvecs + 3 * newton_its + 2 * steps))
        return obs, cost, wall, dict(matvecs=matvecs, newton_its=newton_its,
                                     newton_fail=int(newton_fail), gmres_fail=int(gmres_fail),
                                     converged=newton_fail == 0 and gmres_fail == 0)

    def runner(self, observable: str):
        def f(theta):
            obs, cost, _, _ = self.run(theta)
            return obs[observable], cost

        return f
