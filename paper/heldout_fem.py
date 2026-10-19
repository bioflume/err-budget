"""Held-out problem H1: a variable-coefficient Poisson problem in 3D with a
manufactured solution, solved with scikit-fem (P2 tetrahedra) and an
algebraic-multigrid-preconditioned conjugate-gradient method (pyamg).

    -div(kappa grad u) = f in (0,1)^3,  u = g on the boundary,
    kappa = 1 + x,  u = sin(pi x) sin(pi y) sin(pi z) + x y z.

Observables with exact values:
    mean   = integral of u over the cube                = 8/pi^3 + 1/8
    flux   = integral over the face x = 1 of kappa du/dx = 2 (1/4 - 4/pi)
             (computed from the gradient of the finite-element solution)

Parameters (none of them was used while the program was developed):
    n         cells per side of the hexahedral grid that is split into tetrahedra
                                                               discretisation
    cg_rtol   relative tolerance of the conjugate-gradient solve  nested solver
    amg_theta strength-of-connection threshold of the AMG hierarchy (cost only)

Cost: a deterministic work count, (CG iterations + 1) x nonzeros of the matrix
plus the AMG setup, measured as nonzeros of the hierarchy.
"""
from __future__ import annotations

import math
import time

import numpy as np

EXACT = dict(mean=8.0 / math.pi**3 + 1.0 / 8.0, flux=2.0 * (0.25 - 4.0 / math.pi))
NAMES = ["n", "cg_rtol", "amg_theta"]
LAD = [[3, 4, 6, 8, 12, 16, 24], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [0.0, 0.08, 0.25]]
KINDS = ["discretisation", "nested-solver", "cost-only"]
HS = [[1.0 / n for n in LAD[0]], list(LAD[1]), [1.0 / (1.0 + t) for t in LAD[2]]]


def u_ex(x, y, z):
    return np.sin(np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z) + x * y * z


def f_rhs(x, y, z):
    # -div((1+x) grad u) = -(1+x) lap u - du/dx
    s = np.sin(np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z)
    lap = -3 * np.pi**2 * s
    dudx = np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) * np.sin(np.pi * z) + y * z
    return -(1 + x) * lap - dudx


# At n = 24 the global assembly and the direct solve of the projection of the
# boundary data exceed the 7 GB of the development machine.  From this size on the
# forms are assembled over chunks of elements and the projection is solved by
# conjugate gradients to a relative residual of 1e-14 (verified at n = 8 and 16 to change
# both observables by less than 1e-12 relative; `python heldout_fem.py check`).
LOWMEM_N = 24
_CHUNK = 12000


def _assemble_lowmem(m, a, l):
    import scipy.sparse.linalg as spla
    from skfem import Basis, BilinearForm, ElementTetP2, LinearForm
    from skfem.helpers import dot

    @BilinearForm
    def mass(u, v, w):
        return u * v

    @LinearForm
    def gload(v, w):
        return u_ex(*w.x) * v

    A = b = M = f = None
    nel = m.t.shape[1]
    for s0 in range(0, nel, _CHUNK):
        cb = Basis(m, ElementTetP2(), intorder=4, elements=np.arange(s0, min(s0 + _CHUNK, nel)))
        parts = (a.assemble(cb), l.assemble(cb), mass.assemble(cb), gload.assemble(cb))
        if A is None:
            A, b, M, f = parts
        else:
            A, b, M, f = A + parts[0], b + parts[1], M + parts[2], f + parts[3]
        del cb
    x, info = spla.cg(M.tocsr(), f, rtol=1e-14, atol=0.0, maxiter=5000)
    if info != 0:
        raise RuntimeError("projection of the boundary data did not converge")
    return A.tocsr(), b, x


def run(theta: dict):
    import pyamg
    import scipy.sparse.linalg as spla
    from skfem import (Basis, BilinearForm, ElementTetP2, FacetBasis, Functional,
                       LinearForm, MeshTet, condense)
    from skfem.helpers import dot, grad

    n = int(theta["n"])
    rtol = float(theta["cg_rtol"])
    th = float(theta["amg_theta"])
    t0 = time.perf_counter()
    m = MeshTet.init_tensor(*(np.linspace(0, 1, n + 1),) * 3)
    ib = Basis(m, ElementTetP2(), intorder=4) if n < LOWMEM_N else None

    @BilinearForm
    def a(u, v, w):
        return (1.0 + w.x[0]) * dot(grad(u), grad(v))

    @LinearForm
    def l(v, w):
        return f_rhs(*w.x) * v

    if n < LOWMEM_N:
        A = a.assemble(ib)
        b = l.assemble(ib)
        x = ib.project(lambda X: u_ex(*X))      # boundary values (L2 projection)
    else:
        A, b, x = _assemble_lowmem(m, a, l)
    D = (ib if ib is not None else Basis(m, ElementTetP2(), intorder=1)).get_dofs().all()
    Ai, bi, xi, I = condense(A, b, x=x, D=D)
    ml = pyamg.smoothed_aggregation_solver(Ai.tocsr(), strength=("symmetric", {"theta": th}))
    M = ml.aspreconditioner()
    its = [0]
    sol, info = spla.cg(Ai, bi, rtol=rtol, atol=0.0, M=M, maxiter=2000,
                        callback=lambda _x: its.__setitem__(0, its[0] + 1))
    x[I] = sol
    converged = info == 0
    # observables
    if n < LOWMEM_N:
        mean = float(Functional(lambda w: w["u"]).assemble(ib, u=ib.interpolate(x)))
    else:
        mean, nel = 0.0, m.t.shape[1]
        for s0 in range(0, nel, _CHUNK):
            cb = Basis(m, ElementTetP2(), intorder=4,
                       elements=np.arange(s0, min(s0 + _CHUNK, nel)))
            mean += float(Functional(lambda w: w["u"]).assemble(cb, u=cb.interpolate(x)))
    fb = FacetBasis(m, ElementTetP2(), facets=m.facets_satisfying(lambda X: X[0] > 1 - 1e-12),
                    intorder=4)
    flux = float(Functional(lambda w: (1.0 + w.x[0]) * grad(w["u"])[0]).assemble(
        fb, u=fb.interpolate(x)))
    nnz_h = sum(lv.A.nnz for lv in ml.levels)
    cost = float((its[0] + 1) * Ai.nnz + 2 * nnz_h)
    return (dict(mean=mean, flux=flux), cost, time.perf_counter() - t0,
            dict(cg_its=its[0], converged=bool(converged), dofs=int(A.shape[0])))


def _check(n=16):
    """The low-memory path against the standard one at size n."""
    global LOWMEM_N
    th = dict(n=n, cg_rtol=1e-10, amg_theta=0.0)
    o1 = run(th)[0]
    LOWMEM_N = n
    o2 = run(th)[0]
    LOWMEM_N = 24
    for k in o1:
        print(k, o1[k], o2[k], f"relative difference {abs(o1[k]-o2[k])/abs(o1[k]):.1e}")


if __name__ == "__main__":
    import sys
    if sys.argv[1:2] == ["check"]:
        _check(int(sys.argv[2]) if len(sys.argv) > 2 else 12)
        raise SystemExit
    for n in LAD[0]:
        o, c, w, i = run(dict(n=n, cg_rtol=1e-10, amg_theta=0.0))
        print(n, {k: f"{abs(v-EXACT[k])/abs(EXACT[k]):.2e}" for k, v in o.items()}, f"{c:.3g}",
              f"{w:.1f}s", i)
