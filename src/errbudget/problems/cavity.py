"""L5: a driven cavity with a boundary layer, and observables that are not integrals.

Steady incompressible flow in the unit square in streamfunction form,

    nu Lap(omega) - (psi_y omega_x - psi_x omega_y) + f = 0,   omega = -Lap(psi),

with psi = 0 on all four walls, no slip on three of them, and a moving lid at
y = 1 whose velocity is u = psi_y = g(x).  The manufactured solution

    psi(x, y) = g(x) y^2 q(y),   g(x) = 16 x^2 (1-x)^2,
    q(y) = -delta (1 - exp(-(1-y)/delta)),

satisfies every boundary condition exactly and has a boundary layer of width
delta below the lid, where the velocity rises from about -2 delta to g(x).  The
source f is derived symbolically, so exact values of every observable exist.

Unlike L2-L4 the domain has walls, the error is dominated by a boundary layer,
and none of the three observables is an integral:

    psi_min     the minimum of the streamfunction (an extremum, found on the
                grid and refined by a quadratic through the neighbouring nodes)
    y_min       its height (the location of an extremum, inside the layer)
    omega_lid   the vorticity at the middle of the lid (a point value of a
                second derivative at the wall)

Discretisation: second-order central differences on a uniform N x N grid, with
the wall conditions on psi_n imposed through ghost nodes (the 13-point
biharmonic stencil in disguise).  The nonlinear system is solved by Newton's
method with backtracking from psi = 0; each Newton step by GMRES preconditioned
with an incomplete LU factorisation of the Jacobian.

    N           grid intervals per direction        discretisation
    newton_tol  relative residual reduction         nested solver
    gmres_tol   relative tolerance of each solve    nested solver
    ilu_drop    drop tolerance of the ILU           cost-only

Cost is counted in work units, as for L2: for every Newton step the nonzeros of
the incomplete factors (factorisation) plus, per GMRES iteration, the nonzeros
of the Jacobian and of the factors.  It is deterministic.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _manufactured(nu, delta):
    import sympy as sm
    x, y = sm.symbols("x y")
    g = 16 * x**2 * (1 - x)**2
    q = -delta * (1 - sm.exp(-(1 - y) / delta))
    psi = g * y**2 * q
    lap = lambda w: sm.diff(w, x, 2) + sm.diff(w, y, 2)
    om = -lap(psi)
    f = (sm.diff(psi, y) * sm.diff(om, x) - sm.diff(psi, x) * sm.diff(om, y)) - nu * lap(om)
    L = lambda e: sm.lambdify((x, y), e, "numpy")
    return L(psi), L(om), L(f), sm.lambdify(x, g, "numpy")


class DrivenCavity:
    def __init__(self, nu=1e-2, delta=0.05, max_newton=20):
        self.nu, self.delta, self.max_newton = nu, delta, max_newton
        self.psi_ex, self.om_ex, self.f, self.g = _manufactured(nu, delta)

    # ---------------------------------------------------------------- truth
    def exact_observables(self) -> dict:
        from scipy.optimize import minimize_scalar
        d = self.delta
        m = lambda y: -d * y**2 * (1 - np.exp(-(1 - y) / d))   # psi(0.5, y), g(0.5) = 1
        r = minimize_scalar(m, bounds=(0.5, 1.0), method="bounded",
                            options=dict(xatol=1e-14))
        return {"psi_min": float(r.fun), "y_min": float(r.x),
                "omega_lid": float(self.om_ex(0.5, 1.0))}

    # ------------------------------------------------------------ operators
    def _operators(self, N):
        """Affine map psi_I -> full psi and full omega, and interior stencils."""
        n1 = N + 1
        h = 1.0 / N
        idx = np.arange(n1 * n1).reshape(n1, n1)          # [i, j]: x_i, y_j
        inner = idx[1:-1, 1:-1].ravel()
        nI = len(inner)
        P = sp.csr_matrix((np.ones(nI), (inner, np.arange(nI))), shape=(n1 * n1, nI))
        # 5-point Laplacian rows for every node with all four neighbours present
        rows, cols, vals = [], [], []
        for (di, dj, w) in [(0, 0, -4.0), (1, 0, 1.0), (-1, 0, 1.0), (0, 1, 1.0), (0, -1, 1.0)]:
            ii, jj = np.meshgrid(np.arange(1, N), np.arange(1, N), indexing="ij")
            rows.append(idx[ii, jj].ravel()); cols.append(idx[ii + di, jj + dj].ravel())
            vals.append(np.full(ii.size, w / h**2))
        Lap = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                            shape=(n1 * n1, n1 * n1))
        # omega = -Lap psi in the interior; on the walls from the ghost nodes:
        #   bottom/left/right: omega = -2 psi_1 / h^2;  lid: -(2 psi_{N-1} + 2 h g)/h^2
        r2, c2, v2 = [Lap.tocoo().row], [Lap.tocoo().col], [-Lap.tocoo().data]
        b = np.zeros(n1 * n1)
        x = np.linspace(0, 1, n1)
        for i in range(1, N):
            for (wall, nb) in [((i, 0), (i, 1)), ((0, i), (1, i)), ((N, i), (N - 1, i)),
                               ((i, N), (i, N - 1))]:
                r2.append([idx[wall]]); c2.append([idx[nb]]); v2.append([-2.0 / h**2])
            b[idx[i, N]] = -2.0 * self.g(x[i]) / h
        Om = sp.csr_matrix((np.concatenate(v2), (np.concatenate(r2), np.concatenate(c2))),
                           shape=(n1 * n1, n1 * n1))
        # central first differences at interior nodes
        def d1(di, dj):
            ii, jj = np.meshgrid(np.arange(1, N), np.arange(1, N), indexing="ij")
            r = np.concatenate([idx[ii, jj].ravel()] * 2)
            c = np.concatenate([idx[ii + di, jj + dj].ravel(), idx[ii - di, jj - dj].ravel()])
            v = np.concatenate([np.full(ii.size, 0.5 / h), np.full(ii.size, -0.5 / h)])
            return sp.csr_matrix((v, (r, c)), shape=(n1 * n1, n1 * n1))
        Dx, Dy = d1(1, 0), d1(0, 1)
        R = sp.csr_matrix((np.ones(nI), (np.arange(nI), inner)), shape=(nI, n1 * n1))
        return dict(P=P, Lap=Lap, Om=Om, b=b, Dx=Dx, Dy=Dy, R=R, idx=idx, inner=inner, h=h)

    def run(self, theta: dict):
        N = int(theta["N"])
        ntol = float(theta["newton_tol"])
        gtol = float(theta["gmres_tol"])
        drop = float(theta["ilu_drop"])
        op = self._operators(N)
        P, Lap, Om, b, Dx, Dy, R = (op[k] for k in ("P", "Lap", "Om", "b", "Dx", "Dy", "R"))
        h = op["h"]
        x = np.linspace(0, 1, N + 1)
        X, Y = np.meshgrid(x, x, indexing="ij")
        F = self.f(X, Y).ravel()[op["inner"]]
        A = (Om @ P).tocsr()                    # omega_full = A psi_I + b
        LA = (R @ Lap @ A).tocsr()
        RDxP, RDyP = (R @ Dx @ P).tocsr(), (R @ Dy @ P).tocsr()
        RDxA, RDyA = (R @ Dx @ A).tocsr(), (R @ Dy @ A).tocsr()
        Lb = R @ (Lap @ b)
        Dxb, Dyb = R @ (Dx @ b), R @ (Dy @ b)
        nu = self.nu

        def resid(s):
            px, py = RDxP @ s, RDyP @ s
            ox, oy = RDxA @ s + Dxb, RDyA @ s + Dyb
            return nu * (LA @ s + Lb) - (py * ox - px * oy) + F, (px, py, ox, oy)

        def jac(parts):
            px, py, ox, oy = parts
            D = sp.diags
            return (nu * LA - D(ox) @ RDyP - D(py) @ RDxA + D(oy) @ RDxP + D(px) @ RDyA).tocsc()

        nI = len(op["inner"])
        s = np.zeros(nI)
        r, parts = resid(s)
        r0 = np.linalg.norm(r)
        cost, its, newton, converged = 0.0, 0, 0, False
        gmres_fail = 0
        t0 = time.perf_counter()
        while newton < self.max_newton:
            if np.linalg.norm(r) <= ntol * r0:
                converged = True
                break
            J = jac(parts)
            try:
                ilu = spla.spilu(J, drop_tol=drop, fill_factor=20)
                nnzLU = ilu.L.nnz + ilu.U.nnz
                M = spla.LinearOperator(J.shape, ilu.solve)
            except RuntimeError:
                nnzLU, M = 0, None
            cnt = [0]
            dx, info = spla.gmres(J, -r, rtol=gtol, atol=0.0, M=M, restart=50, maxiter=20,
                                  callback=lambda _r: cnt.__setitem__(0, cnt[0] + 1),
                                  callback_type="pr_norm")
            gmres_fail += info > 0
            its += cnt[0]
            cost += nnzLU + cnt[0] * (J.nnz + nnzLU)
            # backtracking on the residual norm
            lam, nr = 1.0, np.linalg.norm(r)
            while True:
                s_new = s + lam * dx
                r_new, parts_new = resid(s_new)
                cost += J.nnz
                if np.linalg.norm(r_new) < (1 - 1e-4 * lam) * nr or lam < 1e-3:
                    break
                lam *= 0.5
            s, r, parts = s_new, r_new, parts_new
            newton += 1
        else:
            converged = np.linalg.norm(r) <= ntol * r0
        wall = time.perf_counter() - t0
        psi = (P @ s).reshape(N + 1, N + 1)
        om = (A @ s + b).reshape(N + 1, N + 1)
        i, j = np.unravel_index(np.argmin(psi), psi.shape)
        i, j = min(max(i, 1), N - 1), min(max(j, 1), N - 1)

        def vertex(fm, f0, fp):
            den = fm - 2 * f0 + fp
            if den <= 0:
                return 0.0, f0
            t = 0.5 * (fm - fp) / den
            return t, f0 - 0.125 * (fm - fp) ** 2 / den
        tx, _ = vertex(psi[i - 1, j], psi[i, j], psi[i + 1, j])
        ty, _ = vertex(psi[i, j - 1], psi[i, j], psi[i, j + 1])
        dx2 = 0.5 * (psi[i - 1, j] - 2 * psi[i, j] + psi[i + 1, j])
        dy2 = 0.5 * (psi[i, j - 1] - 2 * psi[i, j] + psi[i, j + 1])
        bx = 0.5 * (psi[i + 1, j] - psi[i - 1, j])
        by = 0.5 * (psi[i, j + 1] - psi[i, j - 1])
        pmin = psi[i, j] + bx * tx + dx2 * tx**2 + by * ty + dy2 * ty**2
        obs = {"psi_min": float(pmin), "y_min": float((j + ty) * h),
               "omega_lid": float(om[N // 2, N]) if N % 2 == 0 else
               float(0.5 * (om[N // 2, N] + om[N // 2 + 1, N]))}
        return obs, float(cost), wall, dict(newton=newton, gmres=its,
                                            newton_converged=bool(converged),
                                            gmres_fail=int(gmres_fail),
                                            converged=bool(converged) and gmres_fail == 0)

    def runner(self, observable: str):
        def f(theta):
            obs, cost, _, _ = self.run(theta)
            return obs[observable], cost
        return f
