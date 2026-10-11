"""L3: 2D Burgers, where the linear solve dominates the cost.

    u_t + u (u_x + u_y) = nu (u_xx + u_yy) + f      on [0,2pi)^2, periodic

Every earlier rung was 1D with cheap linear algebra, so the cost half of the
claim was never really tested.  Here each implicit step is a sparse 2D system
that has to be preconditioned, and the preconditioner's own parameter joins the
knobs.

Ground truth by the method of manufactured solutions: we choose u, derive f
symbolically, and the observables follow from the chosen u exactly.  No fine
reference run is used as a stand-in for truth anywhere.

Five knobs, and the fifth is a different species:

    N, steps, newton_tol, gmres_tol     trade accuracy against cost
    ilu_drop                            trades cost against cost

A looser drop tolerance makes the factorisation cheap and the iteration count
high; a tighter one does the reverse, so its cost curve has an interior
minimum.  It has essentially no effect on the answer, because the GMRES
stopping criterion -- not the preconditioner -- sets the accuracy.  Users
routinely tune it as though it controlled accuracy.  A method that claims to
say which knob matters should be able to say that this one does not.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _manufactured(nu):
    """u(x,y,t) and the forcing f that makes it an exact solution."""
    import sympy as sm

    x, y, t = sm.symbols("x y t")
    u = sm.exp(-nu * t) * (sm.sin(x) * sm.cos(y) + sm.Rational(1, 2) * sm.sin(3 * x) * sm.cos(2 * y))
    f = sm.diff(u, t) + u * (sm.diff(u, x) + sm.diff(u, y)) - nu * (sm.diff(u, x, 2) + sm.diff(u, y, 2))
    ux, uy = sm.diff(u, x), sm.diff(u, y)
    return (sm.lambdify((x, y, t), u, "numpy"),
            sm.lambdify((x, y, t), sm.simplify(f), "numpy"),
            sm.lambdify((x, y, t), ux, "numpy"),
            sm.lambdify((x, y, t), uy, "numpy"))


class Burgers2D:
    def __init__(self, nu: float = 0.05, T: float = 1.0, max_newton: int = 20):
        self.nu, self.T, self.max_newton = nu, T, max_newton
        self.u_ex, self.f_ex, self.ux_ex, self.uy_ex = _manufactured(nu)

    # ---------------------------------------------------------------- truth
    def exact_observables(self, Nf: int = 1024) -> dict:
        x = 2 * np.pi * np.arange(Nf) / Nf
        X, Y = np.meshgrid(x, x, indexing="ij")
        u = self.u_ex(X, Y, self.T)
        ux, uy = self.ux_ex(X, Y, self.T), self.uy_ex(X, Y, self.T)
        dA = (2 * np.pi / Nf) ** 2
        return {"energy": 0.5 * float(np.sum(u**2)) * dA,
                "gradient": float(np.sum(ux**2 + uy**2)) * dA}

    # ------------------------------------------------------------ operators
    @staticmethod
    def _ops(N):
        h = 2 * np.pi / N
        o = np.ones(N)
        D1 = sp.diags([o, -o], [1, -1], shape=(N, N), format="lil") / (2 * h)
        D1[0, -1], D1[-1, 0] = -1 / (2 * h), 1 / (2 * h)
        L1 = sp.diags([o, -2 * o, o], [-1, 0, 1], shape=(N, N), format="lil") / h**2
        L1[0, -1], L1[-1, 0] = 1 / h**2, 1 / h**2
        D1, L1 = D1.tocsr(), L1.tocsr()
        I = sp.identity(N, format="csr")
        Dx, Dy = sp.kron(D1, I, "csr"), sp.kron(I, D1, "csr")
        Lap = sp.kron(L1, I, "csr") + sp.kron(I, L1, "csr")
        return Dx, Dy, Lap

    def _F(self, u, Dx, Dy, Lap, fvec):
        return -u * (Dx @ u + Dy @ u) + self.nu * (Lap @ u) + fvec

    def _J(self, u, Dx, Dy, Lap, dt, n):
        Fp = (-sp.diags(Dx @ u + Dy @ u) - sp.diags(u) @ (Dx + Dy)
              + self.nu * Lap)
        return (sp.identity(n, format="csr") - 0.5 * dt * Fp).tocsr()

    # ------------------------------------------------------------------ run
    def run(self, theta: dict):
        N = int(theta["N"])
        steps = int(theta["steps"])
        ntol = float(theta["newton_tol"])
        gtol = float(theta["gmres_tol"])
        drop = float(theta["ilu_drop"])
        # three further knobs, defaulted so earlier studies are unaffected
        restart = int(theta.get("gmres_restart", 50))
        pfreq = int(theta.get("precond_freq", 1))     # refactorise every pfreq steps
        nmax = int(theta.get("newton_max", self.max_newton))
        n, dt = N * N, self.T / steps

        x = 2 * np.pi * np.arange(N) / N
        X, Y = np.meshgrid(x, x, indexing="ij")
        xf, yf = X.ravel(), Y.ravel()
        u = self.u_ex(xf, yf, 0.0)
        Dx, Dy, Lap = self._ops(N)
        nnzA = Dx.nnz + Dy.nnz + Lap.nnz

        cost = 0.0
        gmres_its = newton_its = 0
        newton_fail = gmres_fail = 0
        t0 = time.perf_counter()
        for s in range(steps):
            tn, tnp = s * dt, (s + 1) * dt
            un = u.copy()
            Fn = self._F(un, Dx, Dy, Lap, self.f_ex(xf, yf, tn))
            fnp = self.f_ex(xf, yf, tnp)
            scale = np.linalg.norm(un) + 1e-30

            # preconditioner frozen over pfreq steps, as production codes do
            if s % pfreq == 0:
                J0 = self._J(un, Dx, Dy, Lap, dt, n)
                try:
                    ilu = spla.spilu(J0.tocsc(), drop_tol=drop, fill_factor=20)
                    nnzLU = ilu.L.nnz + ilu.U.nnz
                    M = spla.LinearOperator((n, n), ilu.solve)
                except RuntimeError:
                    ilu, nnzLU, M = None, 0, None
                cost += 4.0 * nnzLU                  # factorisation

            v = un.copy()
            for _ in range(nmax):
                G = v - un - 0.5 * dt * (self._F(v, Dx, Dy, Lap, fnp) + Fn)
                newton_its += 1
                cost += nnzA
                if np.linalg.norm(G) <= ntol * scale:
                    break
                J = self._J(v, Dx, Dy, Lap, dt, n)
                cnt = [0]
                cb = lambda _r: cnt.__setitem__(0, cnt[0] + 1)
                dv, info = spla.gmres(J, -G, rtol=gtol, atol=0.0, M=M,
                                      restart=min(restart, n), maxiter=300,
                                      callback=cb, callback_type="pr_norm")
                gmres_fail += info > 0
                it = max(cnt[0], 1)
                gmres_its += it
                cost += it * (nnzA + 2.0 * nnzLU)    # matvec + preconditioner solve
                v = v + dv
            else:
                G = v - un - 0.5 * dt * (self._F(v, Dx, Dy, Lap, fnp) + Fn)
                newton_fail += np.linalg.norm(G) > ntol * scale
            u = v
        wall = time.perf_counter() - t0

        h = 2 * np.pi / N
        ux, uy = Dx @ u, Dy @ u
        obs = {"energy": 0.5 * float(np.sum(u**2)) * h**2,
               "gradient": float(np.sum(ux**2 + uy**2)) * h**2}
        return obs, float(cost), wall, dict(gmres=gmres_its, newton=newton_its,
                                            nnzLU=nnzLU, newton_fail=int(newton_fail),
                                            gmres_fail=int(gmres_fail),
                                            converged=newton_fail == 0 and gmres_fail == 0)

    def runner(self, observable):
        def f(theta):
            obs, cost, _, _ = self.run(theta)
            return obs[observable], cost
        return f
