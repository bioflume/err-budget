"""Held-out problem H2: the lid-driven cavity at Re = 1000 in OpenFOAM simpleFoam,
against the benchmark of Botella & Peyret (1998).

Unit square, lid velocity 1 at y = 1, nu = 1e-3, uniform N x N grid, SIMPLEC,
second-order upwind-biased convection.  Observables (extrema, not integrals):
    u_min   minimum of u on the vertical centreline x = 1/2     -0.3885698
    v_max   maximum of v on the horizontal centreline y = 1/2    0.3769447
    v_min   minimum of v on the horizontal centreline y = 1/2   -0.5270771
Values on a centreline are the average of the two cell columns (rows) beside
it; the extremum is refined by the parabola through the three samples around it.

Parameters:
    N          cells per side                                  discretisation
    res_ctrl   SIMPLE residualControl for p (U: one tenth)       nested solver
    lin_tol    absolute tolerance of the GAMG pressure solver    nested solver
               (U solver ten times it; relTol 0.1 as in the tutorials)
Cost: cells x outer iterations (deterministic).
A run that does not meet res_ctrl within MAX_ITER iterations is reported as not
converged.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time

import numpy as np

from openfoam_cylinder import ENV, of_cmd, _w, _internal_vectors

REF = dict(u_min=-0.3885698, v_max=0.3769447, v_min=-0.5270771)
NAMES = ["N", "res_ctrl", "lin_tol"]
LAD = [[24, 32, 48, 64, 96, 128, 192], [1e-3, 1e-4, 1e-5, 1e-6, 1e-7], [1e-4, 1e-6, 1e-8]]
KINDS = ["discretisation", "nested-solver", "nested-solver"]
HS = [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2])]
MAX_ITER = 30000
NU = 1e-3

SCHEMES = """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U);
             div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""


def block_mesh(N):
    return f"""convertToMeters 1;
vertices ((0 0 0) (1 0 0) (1 1 0) (0 1 0) (0 0 0.1) (1 0 0.1) (1 1 0.1) (0 1 0.1));
blocks (hex (0 1 2 3 4 5 6 7) ({N} {N} 1) simpleGrading (1 1 1));
edges ();
boundary
(
    lid {{ type wall; faces ((3 7 6 2)); }}
    walls {{ type wall; faces ((0 4 7 3) (2 6 5 1) (1 5 4 0)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
"""


def fv_solution(res, lin):
    return f"""solvers
{{
    p {{ solver GAMG; smoother GaussSeidel; tolerance {lin:g}; relTol 0.1; }}
    U {{ solver smoothSolver; smoother symGaussSeidel; tolerance {10*lin:g}; relTol 0.1; }}
}}
SIMPLE
{{
    nNonOrthogonalCorrectors 0; consistent yes;
    pRefCell 0; pRefValue 0;
    residualControl {{ p {res:g}; U {res/10:g}; }}
}}
relaxationFactors {{ equations {{ U 0.9; }} fields {{ p 1; }} }}
"""


def control_dict(end):
    return f"""application simpleFoam; startFrom latestTime; startTime 0; stopAt endTime;
endTime {end}; deltaT 1; writeControl timeStep; writeInterval {end}; purgeWrite 0;
writeFormat ascii; writePrecision 16; writeCompression off; timeFormat general;
timePrecision 8; runTimeModifiable false;
"""


def _extremum(x, y, kind):
    i = int(np.argmin(y) if kind == "min" else np.argmax(y))
    if 0 < i < len(y) - 1:
        den = y[i - 1] - 2 * y[i] + y[i + 1]
        if (den > 0 and kind == "min") or (den < 0 and kind == "max"):
            return float(y[i] - 0.125 * (y[i - 1] - y[i + 1]) ** 2 / den)
    return float(y[i])


class Cavity:
    def __init__(self, workdir=None):
        self.root = workdir or tempfile.mkdtemp(prefix="errbudget_cavity_")
        self.meshes = {}

    def _sh(self, cmd, cwd):
        r = subprocess.run(of_cmd(cmd), cwd=cwd, env=ENV, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{cmd} failed: {r.stdout[-1500:]} {r.stderr[-1500:]}")
        return r.stdout

    def mesh(self, N):
        if N in self.meshes:
            return self.meshes[N]
        d = os.path.join(self.root, f"mesh_{N}")
        shutil.rmtree(d, ignore_errors=True)
        _w(f"{d}/system/blockMeshDict", "dictionary", "blockMeshDict", block_mesh(N))
        _w(f"{d}/system/controlDict", "dictionary", "controlDict", control_dict(MAX_ITER))
        _w(f"{d}/system/fvSchemes", "dictionary", "fvSchemes", SCHEMES)
        _w(f"{d}/system/fvSolution", "dictionary", "fvSolution", fv_solution(1e-5, 1e-6))
        _w(f"{d}/constant/transportProperties", "dictionary", "transportProperties",
           f"transportModel Newtonian;\nnu {NU};\n")
        _w(f"{d}/constant/turbulenceProperties", "dictionary", "turbulenceProperties",
           "simulationType laminar;\n")
        _w(f"{d}/0/U", "volVectorField", "U", """dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    lid { type fixedValue; value uniform (1 0 0); }
    walls { type noSlip; }
    frontAndBack { type empty; }
}
""")
        _w(f"{d}/0/p", "volScalarField", "p", """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    lid { type zeroGradient; }
    walls { type zeroGradient; }
    frontAndBack { type empty; }
}
""")
        self._sh(["blockMesh"], d)
        self.meshes[N] = d
        return d

    def run(self, theta):
        N = int(theta["N"])
        md = self.mesh(N)
        d = os.path.join(self.root, "run")
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(md, d)
        _w(f"{d}/system/fvSolution", "dictionary", "fvSolution",
           fv_solution(float(theta["res_ctrl"]), float(theta["lin_tol"])))
        t0 = time.perf_counter()
        r = subprocess.run(of_cmd(["simpleFoam"]), cwd=d, env=ENV, capture_output=True, text=True)
        wall = time.perf_counter() - t0
        nan = {k: float("nan") for k in REF}
        if r.returncode != 0:
            return nan, float("inf"), wall, {"unavailable": (r.stdout + r.stderr)[-300:]}
        its = len(re.findall(r"^Time = ", r.stdout, flags=re.M))
        conv = "SIMPLE solution converged" in r.stdout
        times = [t for t in os.listdir(d) if re.fullmatch(r"[\d.eE+-]+", t) and t != "0"]
        U = _internal_vectors(os.path.join(d, max(times, key=float), "U"))
        G = U.reshape(N, N, 3)                      # [j (y), i (x)] for blockMesh order
        h = 1.0 / N
        c = (np.arange(N) + 0.5) * h
        if N % 2 == 0:
            ucl = 0.5 * (G[:, N // 2 - 1, 0] + G[:, N // 2, 0])   # u(x = 1/2, y)
            vcl = 0.5 * (G[N // 2 - 1, :, 1] + G[N // 2, :, 1])   # v(x, y = 1/2)
        else:
            ucl, vcl = G[:, N // 2, 0], G[N // 2, :, 1]
        obs = dict(u_min=_extremum(c, ucl, "min"), v_max=_extremum(c, vcl, "max"),
                   v_min=_extremum(c, vcl, "min"))
        shutil.rmtree(d, ignore_errors=True)
        return obs, float(N * N * its), wall, dict(iters=its, converged=conv, cells=N * N)


if __name__ == "__main__":
    p = Cavity()
    for N in LAD[0][:4]:
        o, c, w, i = p.run(dict(N=N, res_ctrl=1e-7, lin_tol=1e-8))
        print(N, {k: f"{abs(v-REF[k])/abs(REF[k]):.2e}" for k, v in o.items()}, c, f"{w:.1f}s", i,
              flush=True)
