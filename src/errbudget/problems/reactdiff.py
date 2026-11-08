"""L4: two-species reaction-diffusion (Gray-Scott kinetics) in 2D.

    u_t = D_u (u_xx + u_yy) - u v^2 + F (1 - u) + s_u(x,y,t)
    v_t = D_v (v_xx + v_yy) + u v^2 - (F + k) v + s_v(x,y,t)

on [0, 2pi)^2 with periodic boundaries.  This is the canonical model of chemical
pattern formation and of activator-inhibitor systems in mathematical biology.

Why it earns its place alongside the Burgers problems: the solver architecture is
different, not just the equation.  There is no Newton iteration.  Instead the
operators are split -- Strang splitting of reaction against diffusion -- with the
reaction sub-cycled by an explicit integrator and the diffusion solved implicitly.
That gives a parameter set of a different shape:

    N            grid points per direction        discretisation
    steps        outer (splitting) time steps     discretisation
    n_react      reaction sub-steps per half-step discretisation
    gmres_tol    tolerance of the diffusion solve nested solver
    ilu_drop     preconditioner drop tolerance    cost-only

so the accuracy is limited by three interacting discretisation errors -- spatial,
splitting, and reaction sub-cycling -- rather than by one mesh.

Ground truth is exact by the method of manufactured solutions: u and v are chosen
and the source terms s_u, s_v that make them exact solutions are derived
symbolically.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .burgers1d import periodic_differences


def _manufactured(Du, Dv, F, k, rate):
    """Return callables for u, v, their sources, and grad v."""
    import sympy as sm

    x, y, t = sm.symbols("x y t")
    u = sm.Rational(1, 2) + sm.Rational(1, 4) * sm.sin(x) * sm.cos(y) * sm.exp(-t / 2)
    v = sm.Rational(1, 4) + sm.Rational(1, 5) * sm.cos(2 * x) * sm.sin(y) * sm.exp(-t / 2)
    lap = lambda w: sm.diff(w, x, 2) + sm.diff(w, y, 2)
    su = sm.diff(u, t) - Du * lap(u) - rate * (-u * v**2 + F * (1 - u))
    sv = sm.diff(v, t) - Dv * lap(v) - rate * (u * v**2 - (F + k) * v)
    L = lambda e: sm.lambdify((x, y, t), sm.simplify(e), "numpy")
    return L(u), L(v), L(su), L(sv), L(sm.diff(v, x)), L(sm.diff(v, y))


class ReactionDiffusion2D:
    def __init__(self, Du=2e-2, Dv=1e-2, F=0.040, k=0.062, T=1.0, rate=20.0):
        self.Du, self.Dv, self.F, self.k, self.T = Du, Dv, F, k, T
        self.rate = rate
        (self.u_ex, self.v_ex, self.su, self.sv,
         self.vx_ex, self.vy_ex) = _manufactured(Du, Dv, F, k, rate)

    # ---------------------------------------------------------------- truth
    def exact_observables(self, Nf: int = 512) -> dict:
        x = 2 * np.pi * np.arange(Nf) / Nf
        X, Y = np.meshgrid(x, x, indexing="ij")
        v = self.v_ex(X, Y, self.T)
        vx, vy = self.vx_ex(X, Y, self.T), self.vy_ex(X, Y, self.T)
        dA = (2 * np.pi / Nf) ** 2
        return {"inhibitor_mass": float(np.sum(v)) * dA,
                "pattern_sharpness": float(np.sum(vx**2 + vy**2)) * dA}

    # ------------------------------------------------------------ operators
    @staticmethod
    def _laplacian(N):
        _, L1 = periodic_differences(N)
        I = sp.identity(N, format="csr")
        return sp.kron(L1, I, "csr") + sp.kron(I, L1, "csr")

    def _react(self, u, v, xf, yf, t0, dt, n_sub):
        """Explicit RK2 sub-cycling of the reaction (and MMS source) terms."""
        h = dt / n_sub
        t = t0
        for _ in range(n_sub):
            def rhs(uu, vv, tt):
                uv2 = uu * vv * vv
                r = self.rate
                return (r * (-uv2 + self.F * (1 - uu)) + self.su(xf, yf, tt),
                        r * (uv2 - (self.F + self.k) * vv) + self.sv(xf, yf, tt))
            k1u, k1v = rhs(u, v, t)
            k2u, k2v = rhs(u + h * k1u, v + h * k1v, t + h)
            u = u + 0.5 * h * (k1u + k2u)
            v = v + 0.5 * h * (k1v + k2v)
            t += h
        return u, v

    # ------------------------------------------------------------------ run
    def run(self, theta: dict):
        N = int(theta["N"])
        steps = int(theta["steps"])
        n_react = int(theta["n_react"])
        gtol = float(theta["gmres_tol"])
        drop = float(theta["ilu_drop"])
        n, dt = N * N, self.T / steps

        x = 2 * np.pi * np.arange(N) / N
        X, Y = np.meshgrid(x, x, indexing="ij")
        xf, yf = X.ravel(), Y.ravel()
        u = self.u_ex(xf, yf, 0.0)
        v = self.v_ex(xf, yf, 0.0)

        L = self._laplacian(N)
        I = sp.identity(n, format="csr")
        cost, its = 0.0, 0
        gmres_fail = 0
        nnzA = L.nnz
        ops = {}
        for tag, D in (("u", self.Du), ("v", self.Dv)):
            A = (I - 0.5 * dt * D * L).tocsr()
            B = (I + 0.5 * dt * D * L).tocsr()
            try:
                ilu = spla.spilu(A.tocsc(), drop_tol=drop, fill_factor=20)
                nnzLU = ilu.L.nnz + ilu.U.nnz
                M = spla.LinearOperator((n, n), ilu.solve)
            except RuntimeError:
                nnzLU, M = 0, None
            cost += 4.0 * nnzLU                     # one factorisation per species
            ops[tag] = (A, B, M, nnzLU)

        t0 = time.perf_counter()
        for s in range(steps):
            tn = s * dt
            u, v = self._react(u, v, xf, yf, tn, 0.5 * dt, n_react)
            cost += 6.0 * n * n_react
            for tag, w in (("u", u), ("v", v)):
                A, B, M, nnzLU = ops[tag]
                cnt = [0]
                cb = lambda _r: cnt.__setitem__(0, cnt[0] + 1)
                sol, info = spla.gmres(A, B @ w, x0=w, rtol=gtol, atol=0.0, M=M,
                                       restart=min(50, n), maxiter=300,
                                       callback=cb, callback_type="pr_norm")
                gmres_fail += info > 0
                it = max(cnt[0], 1)
                its += it
                cost += it * (nnzA + 2.0 * nnzLU)
                if tag == "u":
                    u = sol
                else:
                    v = sol
            u, v = self._react(u, v, xf, yf, tn + 0.5 * dt, 0.5 * dt, n_react)
            cost += 6.0 * n * n_react
        wall = time.perf_counter() - t0

        h = 2 * np.pi / N
        Vg = v.reshape(N, N)
        vx = (np.roll(Vg, -1, 0) - np.roll(Vg, 1, 0)) / (2 * h)
        vy = (np.roll(Vg, -1, 1) - np.roll(Vg, 1, 1)) / (2 * h)
        obs = {"inhibitor_mass": float(np.sum(Vg)) * h**2,
               "pattern_sharpness": float(np.sum(vx**2 + vy**2)) * h**2}
        return obs, float(cost), wall, dict(gmres=its, gmres_fail=int(gmres_fail),
                                            converged=gmres_fail == 0)
