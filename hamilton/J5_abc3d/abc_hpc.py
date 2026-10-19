"""J5 (held-out): decaying Arnold-Beltrami-Childress flow in a periodic box, 3D,
with OpenFOAM pimpleFoam (current release).  An exact solution of the
incompressible Navier-Stokes equations:

    u(x, t) = (sin z + cos y, sin x + cos z, sin y + cos x) exp(-nu t),
    p = -|u|^2 / 2 (+ const),    in (0, 2 pi)^3, all boundaries cyclic,

because the flow is a Beltrami field (curl u = u) whose nonlinear term is a
gradient.  nu = 0.1, T = 2.  Observables at T (exact values):

    ek      mean kinetic energy  (1/V) int |u|^2/2        = 1.5 exp(-2 nu T)
    eps     dissipation rate     nu (1/V) int |grad u|^2   = 3 nu exp(-2 nu T)
    umax    maximum speed (refined by parabolas through the cell values)
                                                          = sqrt(6) exp(-nu T)

Parameters (fixed before any run; none was used in developing the program):
    N         cells per side                          16 24 32 48 64 96     discretisation
    dt        time step                               0.04 0.02 0.01 0.005  discretisation
    p_tol     GAMG absolute tolerance (U: ten times)  1e-3 1e-5 1e-7        nested solver
    n_outer   PIMPLE outer correctors per step        1 2 3                 nested solver
Cost: core-seconds (ExecutionTime x cores).
    python -m farm abc_hpc array K NK ;  python analyse_j5.py
"""
from __future__ import annotations

import itertools
import json
import math
import os
import re
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
import ofmodern as OF
from openfoam_cylinder import _w, _internal_vectors

NU, T_END, L = 0.1, 2.0, 2 * math.pi
EXACT = dict(ek=1.5 * math.exp(-2 * NU * T_END), eps=3 * NU * math.exp(-2 * NU * T_END),
             umax=math.sqrt(6.0) * math.exp(-NU * T_END))
NAMES = ["N", "dt", "p_tol", "n_outer"]
LAD = [[16, 24, 32, 48, 64, 96], [0.04, 0.02, 0.01, 0.005], [1e-3, 1e-5, 1e-7], [1, 2, 3]]
KINDS = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
HS = [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), [1.0 / n for n in LAD[3]]]
NPROCS = {16: 1, 24: 1, 32: 1, 48: 2, 64: 4, 96: 8}
SMOKE = os.environ.get("KNOBS_SMOKE") == "1"

SCHEMES = """ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss linear;
             div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""


def solution(p_tol, n_outer):
    return f"""solvers
{{
    p {{ solver GAMG; smoother GaussSeidel; tolerance {p_tol:g}; relTol 0.01; }}
    pFinal {{ $p; relTol 0; }}
    "(U|UFinal)" {{ solver smoothSolver; smoother symGaussSeidel; tolerance {10*p_tol:g}; relTol 0; }}
}}
PIMPLE {{ momentumPredictor yes; nOuterCorrectors {n_outer}; nCorrectors 2;
          nNonOrthogonalCorrectors 0; pRefCell 0; pRefValue 0; }}
relaxationFactors {{ fields {{ p 1; pFinal 1; }} equations {{ U 1; UFinal 1; }} }}
"""


def block_mesh(N):
    return f"""convertToMeters 1;
vertices ((0 0 0) ({L} 0 0) ({L} {L} 0) (0 {L} 0) (0 0 {L}) ({L} 0 {L}) ({L} {L} {L}) (0 {L} {L}));
blocks (hex (0 1 2 3 4 5 6 7) ({N} {N} {N}) simpleGrading (1 1 1));
edges ();
boundary
(
    left   {{ type cyclic; neighbourPatch right;  faces ((0 4 7 3)); }}
    right  {{ type cyclic; neighbourPatch left;   faces ((1 2 6 5)); }}
    bottom {{ type cyclic; neighbourPatch top;    faces ((0 1 5 4)); }}
    top    {{ type cyclic; neighbourPatch bottom; faces ((3 7 6 2)); }}
    front  {{ type cyclic; neighbourPatch back;   faces ((0 3 2 1)); }}
    back   {{ type cyclic; neighbourPatch front;  faces ((4 5 6 7)); }}
);
"""


def centres(N):
    h = L / N
    c = (np.arange(N) + 0.5) * h
    Z, Y, X = np.meshgrid(c, c, c, indexing="ij")       # x fastest (blockMesh order)
    return X.ravel(), Y.ravel(), Z.ravel()


def u0(x, y, z):
    return np.stack([np.sin(z) + np.cos(y), np.sin(x) + np.cos(z), np.sin(y) + np.cos(x)], axis=1)


def field_text(cls, dims, values, boundary):
    if cls == "volVectorField":
        body = "\n".join(f"({a:.15g} {b:.15g} {c:.15g})" for a, b, c in values)
        tp = "List<vector>"
    else:
        body = "\n".join(f"{a:.15g}" for a in values)
        tp = "List<scalar>"
    return (f"dimensions {dims};\ninternalField nonuniform {tp}\n{len(values)}\n(\n{body}\n);\n"
            f"boundaryField\n{{\n{boundary}}}\n")


CYC = "".join(f"    {p} {{ type cyclic; }}\n" for p in ("left", "right", "bottom", "top", "front", "back"))


def TASKS():
    T = []
    for idx in itertools.product(*[range(len(l)) for l in LAD]):
        th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
        if SMOKE and not (th["N"] == 16 and th["dt"] == 0.04 and th["p_tol"] == 1e-5):
            continue
        T.append(dict(id="abc_" + "_".join(map(str, idx)), idx=list(idx), theta=th,
                      nprocs=NPROCS[th["N"]]))
    return T


def observables(U, N):
    h = L / N
    ek = 0.5 * float(np.mean(np.sum(U**2, axis=1)))
    G = U.reshape(N, N, N, 3)                         # [z, y, x, comp]
    g2 = 0.0
    for ax in range(3):                               # periodic central differences
        d = (np.roll(G, -1, axis=ax) - np.roll(G, 1, axis=ax)) / (2 * h)
        g2 += np.sum(d**2, axis=3)
    eps = NU * float(np.mean(g2))
    s = np.sqrt(np.sum(G**2, axis=3))
    k = np.unravel_index(np.argmax(s), s.shape)
    val = s[k]
    for ax in range(3):                               # parabola along each axis
        km, kp = list(k), list(k)
        km[ax] = (k[ax] - 1) % N; kp[ax] = (k[ax] + 1) % N
        fm, fp = s[tuple(km)], s[tuple(kp)]
        den = fm - 2 * s[k] + fp
        if den < 0:
            val += -0.125 * (fm - fp) ** 2 / den
    return dict(ek=ek, eps=eps, umax=float(val))


def run(task, out):
    th = task["theta"]
    N, dt = int(th["N"]), float(th["dt"])
    T = 0.4 if SMOKE else T_END
    case = os.path.join(out, "case")
    shutil.rmtree(case, ignore_errors=True)
    _w(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", block_mesh(N))
    _w(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", SCHEMES)
    _w(f"{case}/system/fvSolution", "dictionary", "fvSolution", solution(th["p_tol"], th["n_outer"]))
    _w(f"{case}/system/controlDict", "dictionary", "controlDict",
       OF.control_dict(0, T, dt, T))
    for name, body in [("transportProperties", f"transportModel Newtonian;\nnu {NU};\n"),
                       ("turbulenceProperties", "simulationType laminar;\n"),
                       ("momentumTransport", "simulationType laminar;\n")]:
        _w(f"{case}/constant/{name}", "dictionary", name, body)
    OF._sh(["blockMesh"], case)
    x, y, z = centres(N)
    U = u0(x, y, z)
    _w(f"{case}/0/U", "volVectorField", "U", field_text("volVectorField", "[0 1 -1 0 0 0 0]", U, CYC))
    _w(f"{case}/0/p", "volScalarField", "p",
       field_text("volScalarField", "[0 2 -2 0 0 0 0]", -0.5 * np.sum(U**2, axis=1), CYC))
    log = os.path.join(case, "log")
    n = task["nprocs"] if not SMOKE else 1
    failed = None
    try:
        OF.run("pimpleFoam", case, nprocs=n, log=log)
        if n > 1:
            OF._sh(["reconstructPar", "-latestTime"], case)
    except RuntimeError as exc:
        failed = str(exc)[-600:]
    steps = OF.log_steps(log)
    res = dict(theta=th, nprocs=n, failed=failed, steps=len(steps))
    last = OF.latest_time(case)
    if failed is None and last is not None and abs(float(last) - T) < 1e-9:
        Uf = _internal_vectors(os.path.join(case, last, "U"))
        res["obs"] = observables(Uf, N)
        res["cost"] = float(steps[-1][1] * n) if steps else None
        res["err"] = {k: abs(v - EXACT[k] * (math.exp(-2 * NU * (T - T_END)) if k != "umax"
                                             else math.exp(-NU * (T - T_END)))) / EXACT[k]
                      for k, v in res["obs"].items()}
    elif len(steps) == 0 and failed:
        raise RuntimeError(f"pimpleFoam did not start:\n{failed}")
    json.dump(res, open(os.path.join(out, "result.json"), "w"), indent=1)
    shutil.rmtree(case, ignore_errors=True)
