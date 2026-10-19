"""DFG benchmark "flow around a cylinder", case 2D-2 (periodic, Re = 100), OpenFOAM.

Same geometry and mesh family as 2D-1 (openfoam_cylinder.py), with a steady
parabolic inflow of maximum velocity 1.5 (mean 1.0), so Re = 100 and the flow
sheds vortices with a period of about 0.33 s.  The observables are the ones of
the benchmark: the maximum drag and lift coefficients over a period and the
Strouhal number St = f D / U_mean.

The solver is OpenFOAM's pimpleFoam (PIMPLE, second-order backward time
stepping, fixed time step), driven only through its input files:

    n_arc       cells per quarter of the cylinder (whole mesh scales)   discretisation
    dt          time step                                               discretisation
    outer_tol   PIMPLE outerCorrectorResidualControl for p and U; the
                outer iteration within a step stops there (at most 50)  nested solver
    p_tol       absolute tolerance of the GAMG pressure solver (U solver:
                ten times it)                                           nested solver

Every run starts from a periodic state of its own mesh (one spin-up per mesh
level, 0 -> T_SPIN with the finest settings), relaxes for T_RELAX (2 s) so that the
start state is forgotten, and is then sampled every SAMPLE s for T_MEAS (2 s,
about six shedding periods).  Spin-up is 20 s from rest.  The
Ubuntu package of OpenFOAM v1912 cannot construct function objects, so the
fields are written at the sampling times and the force coefficients are
computed from them as the forces function object does (see
openfoam_cylinder.py).  Cost is the CPU time OpenFOAM reports for the relax
and measure phases (the spin-up is shared by all runs on a mesh and excluded).
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time

import numpy as np

from openfoam_cylinder import (ENV, of_cmd, H, D, DEPTH, NU, _w, block_mesh_dict, _patch_centres,
                               _wall_geometry, _internal_vectors, _internal_scalars)

UM2 = 1.5                      # maximum inflow velocity of 2D-2
UBAR2 = 2.0 * UM2 / 3.0        # = 1.0
T_SPIN = float(os.environ.get("KNOBS_DFG2_TSPIN", 20.0))
T_RELAX = float(os.environ.get("KNOBS_DFG2_TRELAX", 2.0))
T_MEAS = float(os.environ.get("KNOBS_DFG2_TMEAS", 2.0))
RELAX_CHECKPOINT = 0.5          # checkpoint interval of the relaxation phase
SAMPLE = 0.004                 # a multiple of every time step on the ladder
# Published values for orientation only (Schaefer & Turek 1996 intervals; the
# FEATFLOW level-6 run): grading uses a converged OpenFOAM reference whose own
# uncertainty is measured (exp56_dfg2d2.py).
PUBLISHED = dict(cd_max=(3.22, 3.24), cl_max=(0.99, 1.01), st=(0.295, 0.305))
# FEATFLOW, level 6 (667,264 unknowns), time step 1/200: the finest values on the
# benchmark page (featflow.de, "DFG benchmark 2D-2"); level 5 at 1/400 gives
# c_D,max 3.2122, so these are themselves uncertain by a few parts in 1e3
FEATFLOW = dict(cd_max=3.2200, cl_max=0.9859, st=0.30188)

FV_SCHEMES = """ddtSchemes { default backward; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""


def fv_solution(outer_tol, p_tol, n_outer=50):
    return f"""solvers
{{
    p {{ solver GAMG; smoother GaussSeidel; tolerance {p_tol:g}; relTol 0.05; }}
    pFinal {{ $p; relTol 0; }}
    "(U|UFinal)" {{ solver smoothSolver; smoother symGaussSeidel; tolerance {10*p_tol:g}; relTol 0; }}
}}
PIMPLE
{{
    momentumPredictor yes;
    nOuterCorrectors {n_outer};
    nCorrectors 2;
    nNonOrthogonalCorrectors 0;
    residualControl
    {{
        p {{ tolerance {outer_tol:g}; relTol 0; }}
        U {{ tolerance {outer_tol:g}; relTol 0; }}
    }}
}}
relaxationFactors
{{
    fields {{ p 1; pFinal 1; }}
    equations {{ U 1; UFinal 1; }}
}}
"""


def control_dict(start, end, dt, write_every=None):
    if write_every is None:
        wc = f"writeControl runTime; writeInterval {end - start:.10g};"
    else:
        wc = f"writeControl runTime; writeInterval {write_every:.10g};"
    return f"""application pimpleFoam;
startFrom latestTime; startTime {start:.10g}; stopAt endTime; endTime {end:.10g};
deltaT {dt:.10g}; {wc} purgeWrite 0;
writeFormat ascii; writePrecision 16; writeCompression off;
timeFormat general; timePrecision 10; runTimeModifiable false; adjustTimeStep no;
"""


def _coeffs(tdir, geo):
    U = _internal_vectors(os.path.join(tdir, "U"))[geo["cells"]]
    p = _internal_scalars(os.path.join(tdir, "p"))[geo["cells"]]
    sf, cf, cc = geo["sf"], geo["cf"], geo["cc"]
    mag = np.linalg.norm(sf, axis=1)
    nhat = sf / mag[:, None]
    d = np.abs(np.einsum("ij,ij->i", cc - cf, nhat))
    ut = U - np.einsum("ij,ij->i", U, nhat)[:, None] * nhat
    F = (p[:, None] * sf).sum(axis=0) + (NU * mag[:, None] * ut / d[:, None]).sum(axis=0)
    q = 0.5 * UBAR2**2 * D * DEPTH
    return float(F[0] / q), float(F[1] / q)


def _refined_max(t, y, i):
    """Value at sample i refined by the parabola through its neighbours."""
    if i <= 0 or i >= len(y) - 1:
        return float(y[i])
    den = y[i - 1] - 2 * y[i] + y[i + 1]
    return float(y[i] - 0.125 * (y[i - 1] - y[i + 1]) ** 2 / den) if den < 0 else float(y[i])


def _crossings(t, y):
    """Times (and sample indices) of the upward zero crossings of y - mean(y)."""
    z = y - y.mean()
    out = []
    for i in range(len(z) - 1):
        if z[i] < 0 <= z[i + 1]:
            out.append((t[i] - z[i] * (t[i + 1] - t[i]) / (z[i + 1] - z[i]), i))
    return out


def _peak(t, y, ref):
    """Mean over complete shedding periods (between successive upward crossings
    of the lift `ref`) of the refined maximum of y within each period.  Small
    sample-to-sample wiggles cannot create spurious maxima this way."""
    cr = _crossings(t, ref)
    vals = []
    for (_, i0), (_, i1) in zip(cr[:-1], cr[1:]):
        j = i0 + int(np.argmax(y[i0:i1 + 1]))
        vals.append(_refined_max(t, y, j))
    return (float(np.mean(vals)) if vals else float("nan")), len(vals)


def _frequency(t, y):
    """Frequency from the upward zero crossings of y - mean(y)."""
    z = y - y.mean()
    tc = [t[i] - z[i] * (t[i + 1] - t[i]) / (z[i + 1] - z[i])
          for i in range(len(z) - 1) if z[i] < 0 <= z[i + 1]]
    if len(tc) < 2:
        return float("nan")
    return float((len(tc) - 1) / (tc[-1] - tc[0]))


def _times(d, after):
    """Time directories of case d later than `after`, sorted."""
    return sorted((float(x), x) for x in os.listdir(d)
                  if re.fullmatch(r"[\d.eE+-]+", x) and float(x) > after + 1e-9)


def _log_upto(log, upto):
    """CPU time and PIMPLE iteration numbers of a pimpleFoam log, counting only the
    time steps up to `upto` (all of them if None)."""
    steps, cur, its = [], None, []
    for line in open(log, errors="replace"):
        if line.startswith("Time = "):
            cur, its = float(line.split()[2]), []
        elif line.startswith("PIMPLE: iteration"):
            its.append(int(line.split()[2]))
        elif line.startswith("ExecutionTime =") and cur is not None:
            steps.append((cur, float(line.split()[2]), its))
    keep = [x for x in steps if upto is None or x[0] <= upto + 1e-9]
    return (keep[-1][1] if keep else 0.0), [i for x in keep for i in x[2]]


class DFG2D2:
    def __init__(self, workdir=None, spin=None):
        self.root = workdir or tempfile.mkdtemp(prefix="knobs_dfg2_")
        os.makedirs(self.root, exist_ok=True)
        self.spin = spin or dict(dt=0.0005, outer_tol=1e-5, p_tol=1e-8)
        self.meshes = {}

    def _sh(self, cmd, cwd):
        r = subprocess.run(of_cmd(cmd), cwd=cwd, env=ENV, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd)} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
        return r.stdout

    def spin_dt(self, na):
        """Spin-up time step: halved on the meshes finer than n_arc = 32, which keeps
        the maximum Courant number below one (at Courant number 1.3 the spin-up of
        n_arc = 48 diverges)."""
        return self.spin["dt"] * (0.5 if na > 32 else 1.0)

    def spin_outer(self, na):
        """Outer iterations of the spin-up.  On the finest mesh (n_arc = 64) the
        impulsive start from rest makes the unrelaxed outer iteration diverge in the
        first time step, so that spin-up runs in PISO mode (one outer iteration)."""
        return 1 if na > 48 else 50

    def mesh(self, na):
        """Mesh, inflow and a periodic start state for one level (cached on disk)."""
        if na in self.meshes:
            return self.meshes[na]
        d = os.path.join(self.root, f"mesh_{na}")
        done = os.path.join(d, "SPUN_UP")
        setup = os.path.join(d, "SET_UP")
        if not os.path.exists(done) and not os.path.exists(setup):
            shutil.rmtree(d, ignore_errors=True)
            _w(f"{d}/system/blockMeshDict", "dictionary", "blockMeshDict", block_mesh_dict(na))
            _w(f"{d}/system/controlDict", "dictionary", "controlDict",
               control_dict(0.0, T_SPIN, self.spin_dt(na)))
            _w(f"{d}/system/fvSchemes", "dictionary", "fvSchemes", FV_SCHEMES)
            _w(f"{d}/system/fvSolution", "dictionary", "fvSolution",
               fv_solution(self.spin["outer_tol"], self.spin["p_tol"], self.spin_outer(na)))
            _w(f"{d}/constant/transportProperties", "dictionary", "transportProperties",
               f"transportModel Newtonian;\nnu {NU};\n")
            _w(f"{d}/constant/turbulenceProperties", "dictionary", "turbulenceProperties",
               "simulationType laminar;\n")
            self._sh(["blockMesh"], d)
            _, cy = _patch_centres(d, "inlet")
            prof = 4.0 * UM2 * cy * (H - cy) / H**2
            inlet = "nonuniform List<vector>\n%d\n(\n%s\n)" % (
                len(prof), "\n".join(f"({u:.15g} 0 0)" for u in prof))
            _w(f"{d}/0/U", "volVectorField", "U", f"""dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
    inlet {{ type fixedValue; value {inlet}; }}
    outlet {{ type zeroGradient; }}
    walls {{ type noSlip; }}
    cylinder {{ type noSlip; }}
    frontAndBack {{ type empty; }}
}}
""")
            _w(f"{d}/0/p", "volScalarField", "p", """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type zeroGradient; }
    outlet { type fixedValue; value uniform 0; }
    walls { type zeroGradient; }
    cylinder { type zeroGradient; }
    frontAndBack { type empty; }
}
""")
            open(setup, "w").write("mesh and initial fields written\n")
        if not os.path.exists(done):
            # The spin-up writes a checkpoint every second of flow and resumes from the
            # last complete one, so an interrupted spin-up loses at most two seconds.
            # The newest checkpoint may have been cut off mid-write and is discarded.
            times = sorted((float(x), x) for x in os.listdir(d)
                           if re.fullmatch(r"[\d.eE+-]+", x) and float(x) > 0)
            if times:
                shutil.rmtree(os.path.join(d, times[-1][1]))
            _w(f"{d}/system/controlDict", "dictionary", "controlDict",
               control_dict(0.0, T_SPIN, self.spin_dt(na), 1.0))
            _w(f"{d}/system/fvSolution", "dictionary", "fvSolution",
               fv_solution(self.spin["outer_tol"], self.spin["p_tol"], self.spin_outer(na)))
            t0 = time.time()
            self._sh(["pimpleFoam"], d)
            for tv, x in [(float(x), x) for x in os.listdir(d)
                          if re.fullmatch(r"[\d.eE+-]+", x)]:
                if 0 < tv < T_SPIN - 1e-9:
                    shutil.rmtree(os.path.join(d, x))
            open(done, "w").write(f"spin-up 0 -> {T_SPIN} s, last {time.time()-t0:.0f} s wall\n")
        ncell = int(re.search(r"cells:\s+(\d+)", self._sh(["checkMesh"], d)).group(1))
        geo = _wall_geometry(d, "cylinder")
        self.meshes[na] = (d, ncell, geo)
        return self.meshes[na]

    def run(self, theta: dict):
        """theta: n_arc, dt, outer_tol, p_tol -> (obs, cost, wall, info)."""
        md, ncell, geo = self.mesh(int(theta["n_arc"]))
        dt, otol, ptol = float(theta["dt"]), float(theta["outer_tol"]), float(theta["p_tol"])
        # The run directory is named after the setting, and the run checkpoints as it
        # goes (every RELAX_CHECKPOINT of flow while relaxing, every sample while
        # measuring), so a run killed by a machine restart resumes where it stopped.
        # Cost and outer iterations are read from the solver logs up to the checkpoint
        # that is resumed from, so they are those of an uninterrupted run.
        tag = "%d_%g_%g_%g" % (int(theta["n_arc"]), dt, otol, ptol)
        d = os.path.join(self.root, f"run_{tag}")
        state = os.path.join(d, "STATE.json")
        t_s = T_SPIN
        if os.path.exists(state):
            st = json.load(open(state))
        else:
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d)
            for sub in ("constant", "system"):
                shutil.copytree(os.path.join(md, sub), os.path.join(d, sub))
            shutil.copytree(os.path.join(md, f"{t_s:g}"), os.path.join(d, f"{t_s:g}"))
            _w(f"{d}/system/fvSolution", "dictionary", "fvSolution", fv_solution(otol, ptol))
            st = dict(theta=theta, phase=0, seg=0, cpu=0.0, outer=[], wall=0.0)
            json.dump(st, open(state, "w"))
        nan = {"cd_max": float("nan"), "cl_max": float("nan"), "st": float("nan")}
        wall0 = time.perf_counter()
        wall = lambda: st["wall"] + time.perf_counter() - wall0
        phases = [(t_s, t_s + T_RELAX, RELAX_CHECKPOINT),
                  (t_s + T_RELAX, t_s + T_RELAX + T_MEAS, SAMPLE)]
        while st["phase"] < len(phases):
            a, b, every = phases[st["phase"]]
            log = os.path.join(d, "log.%d" % st["seg"])
            if os.path.exists(log):          # this segment was interrupted
                ts = _times(d, a)
                if ts:                       # the newest checkpoint may be incomplete
                    shutil.rmtree(os.path.join(d, ts[-1][1]))
                ts = _times(d, a)
                r = ts[-1][0] if ts else a
                c, o = _log_upto(log, r)
                st["cpu"] += c; st["outer"] += o; st["seg"] += 1
                json.dump(st, open(state, "w"))
                log = os.path.join(d, "log.%d" % st["seg"])
            _w(f"{d}/system/controlDict", "dictionary", "controlDict",
               control_dict(a, b, dt, every))
            with open(log, "w") as fo:
                r = subprocess.run(of_cmd(["pimpleFoam"]), cwd=d, env=ENV, stdout=fo,
                                   stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                last = [l.strip() for l in open(log, errors="replace") if l.startswith("Time =")][-1:]
                err = r.stderr.strip().splitlines()
                shutil.rmtree(d, ignore_errors=True)
                return nan, float("inf"), wall(), {"unavailable": " | ".join(last + err[:4])[:600]}
            c, o = _log_upto(log, None)
            st["cpu"] += c; st["outer"] += o; st["seg"] += 1; st["phase"] += 1
            if st["phase"] == 1:             # drop the relaxation checkpoints
                for tv, x in _times(d, a):
                    if tv < b - 1e-9:
                        shutil.rmtree(os.path.join(d, x))
            st["wall"] = wall(); wall0 = time.perf_counter()
            json.dump(st, open(state, "w"))
        cpu, outer = st["cpu"], st["outer"]
        # force history from the written fields of the measurement window
        times = sorted((float(x), x) for x in os.listdir(d)
                       if re.fullmatch(r"[\d.eE+-]+", x) and float(x) > t_s + T_RELAX - 1e-9)
        t, cd, cl = [], [], []
        for tv, name in times:
            try:
                a, b = _coeffs(os.path.join(d, name), geo)
            except Exception:
                continue
            t.append(tv); cd.append(a); cl.append(b)
        shutil.rmtree(d, ignore_errors=True)
        t, cd, cl = map(np.array, (t, cd, cl))
        if len(t) < 20 or not (np.all(np.isfinite(cd)) and np.all(np.isfinite(cl))):
            return nan, float("inf"), wall(), {"unavailable": "no history"}
        cdm, ncd = _peak(t, cd, cl)
        clm, ncl = _peak(t, cl, cl)
        f = _frequency(t, cl)
        obs = {"cd_max": cdm, "cl_max": clm, "st": f * D / UBAR2}
        info = dict(cells=ncell, samples=len(t), peaks=(ncd, ncl),
                    mean_outer=float(np.mean(outer)) if outer else float("nan"),
                    max_outer=max(outer) if outer else -1,
                    cd_series=cd.tolist(), cl_series=cl.tolist())
        if not all(np.isfinite(list(obs.values()))):
            return nan, float("inf"), wall(), {"unavailable": "no shedding"}
        return obs, cpu, wall(), info
