"""J2: DFG benchmark 2D-2 refined into the asymptotic range, with the averaging
window and the relaxation time as parameters (referee point 8 of revision 3).

Grid of base runs (180 + 3 reference runs):
    n_arc      24, 32, 48, 64, 96                discretisation
    dt         2e-3, 1e-3, 5e-4, 2.5e-4          discretisation
    outer_tol  1e-2, 1e-3, 1e-4                  nested solver (PIMPLE residualControl)
    p_tol      1e-4, 1e-6, 1e-8                  nested solver (GAMG; U ten times it)
Each base run starts from a periodic state of its mesh and integrates T_TOT = 8
time units with OpenFOAM's forceCoeffs function object writing c_D and c_L every
time step.  The relaxation time R in {1, 2, 4} and the averaging window W in
{1, 2, 4} (R + W <= 8) are two further parameters derived from the same run:
the observables are evaluated on [t0 + R, t0 + R + W] and the cost is the CPU
time up to t0 + R + W, which is what a run of that length would cost.  This
gives a 6-parameter tensor of 180 x 9 settings at the cost of 180 runs.

Reference runs: n_arc 192 and 128 at dt = 1.25e-4, and n_arc 128 at 2.5e-4,
with outer_tol 1e-5 and p_tol 1e-8.

Spin-up: n_arc 24 from rest for 20 time units; every finer mesh is initialised
by mapFields from the next coarser spun-up mesh and run for 5 time units.  Spin-ups
run in PISO mode at a time step with maximum Courant number below about 0.9.

Revision 4.  The spin-up time steps of revision 3 of this script (1e-3 at n_arc 24
and 32, halved every two levels) give maximum Courant numbers of about 1.3-1.7 on
this mesh family (the article's 2D-2 driver measured 1.3 at n_arc 48 and 5e-4, where
its spin-up diverged, i.e. Co_max ~ 54 n_arc dt).  They are now the article's values
(Co_max 0.65-0.87); a spin-up that fails is restarted from its initial state at half
the time step, up to twice; the largest Courant number of every spin-up and every run
is recorded; and one mesh level whose spin-up fails no longer stops the others.

November 2026 (rerun on Hamilton).  mapFields -consistent maps the values of every
patch, the fixedValue inlet included, so every mapped mesh started from the coarser
mesh's inlet profile, interpolated, instead of the analytic one (flow rate off by up to
0.4%; the spin-up and runs of revision 4 on Hamilton were affected).  spin() now writes
the analytic inlet back after mapFields and checks it (_restore_inlet), records it in
INLET_ANALYTIC, and spin() and run() refuse a mapped mesh without that record.  SPUN_UP
records the mesh actually mapped from (the MAPPED file), and result.json records
inlet = "analytic".  Measured on this mesh family, Co_max = 27.8 n_arc dt (not 54), so
the spin-up steps below give Co_max 0.33-0.45.  The smoke test spins up its finer mesh
in parallel, as production does.
    python dfg2d2_hpc.py spin [n_arc ...]      # sequential, coarse to fine
    python -m farm dfg2d2_hpc array K NK         # the runs (see submit.sh)
"""
from __future__ import annotations

import glob
import json
import os
import re
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
import ofmodern as OF
from openfoam_cylinder import block_mesh_dict, _patch_centres, H, NU, _w
from openfoam_dfg2d2 import FV_SCHEMES

UM2 = 1.5
NA = [24, 32, 48, 64, 96]
DT = [2e-3, 1e-3, 5e-4, 2.5e-4]
OUTER = [1e-2, 1e-3, 1e-4]
PTOL = [1e-4, 1e-6, 1e-8]
REFS = [dict(n_arc=192, dt=1.25e-4, outer_tol=1e-5, p_tol=1e-8),
        dict(n_arc=128, dt=1.25e-4, outer_tol=1e-5, p_tol=1e-8),
        dict(n_arc=128, dt=2.5e-4, outer_tol=1e-5, p_tol=1e-8)]
NPROCS = {24: 1, 32: 1, 48: 2, 64: 4, 96: 8, 128: 16, 192: 32}
# Co_max ~ 54 n_arc dt (see the docstring): 0.65 at n_arc 24, 48, 96, 192 and 0.87 at 32, 64, 128
SPIN_DT = {24: 5e-4, 32: 5e-4, 48: 2.5e-4, 64: 2.5e-4, 96: 1.25e-4, 128: 1.25e-4, 192: 6.25e-5}
SPIN_RETRIES = 2                 # halvings of the spin-up time step after a failure
T_SPIN0, T_SPIN_MAPPED, T_TOT = 20.0, 5.0, 8.0
SMOKE = os.environ.get("ERRBUDGET_SMOKE") == "1"
if SMOKE:            # a few minutes on the test partition
    NA, DT, OUTER, PTOL = [12, 16], [2e-3], [1e-3], [1e-6]
    REFS = [dict(n_arc=16, dt=2e-3, outer_tol=1e-4, p_tol=1e-8)]
    T_SPIN0, T_SPIN_MAPPED, T_TOT = 1.0, 0.5, 0.8
    NPROCS.update({12: 1, 16: 2})
    SPIN_DT.update({12: 1e-3, 16: 5e-4})
ROOT = os.path.join(os.environ.get("WORK", "."), "dfg2d2_cases")
# relaxation time and averaging window, derived from each run (1, 2, 4 time units)
RELAX = [T_TOT / 8, T_TOT / 4, T_TOT / 2]
WIN = [T_TOT / 8, T_TOT / 4, T_TOT / 2]


def fv_solution(outer_tol, p_tol, n_outer):
    """PIMPLE settings; both spellings of the outer residual control are given
    (ESI releases read residualControl, others outerCorrectorResidualControl;
    an unknown entry is ignored)."""
    rc = f"p {{ tolerance {outer_tol:g}; relTol 0; }} U {{ tolerance {outer_tol:g}; relTol 0; }}"
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
    residualControl {{ {rc} }}
    outerCorrectorResidualControl {{ {rc} }}
}}
relaxationFactors {{ fields {{ p 1; pFinal 1; }} equations {{ U 1; UFinal 1; }} }}
"""


def mesh_dir(na):
    return os.path.join(ROOT, f"mesh_{na}")


def make_mesh(na):
    d = mesh_dir(na)
    if os.path.exists(os.path.join(d, "MESHED")):
        return d
    shutil.rmtree(d, ignore_errors=True)
    _w(f"{d}/system/blockMeshDict", "dictionary", "blockMeshDict", block_mesh_dict(na))
    _w(f"{d}/system/controlDict", "dictionary", "controlDict", OF.control_dict(0, 1, 1e-3, 1))
    _w(f"{d}/system/fvSchemes", "dictionary", "fvSchemes", FV_SCHEMES)
    _w(f"{d}/system/fvSolution", "dictionary", "fvSolution", fv_solution(1e-4, 1e-8, 1))
    _w(f"{d}/constant/transportProperties", "dictionary", "transportProperties",
       f"transportModel Newtonian;\nnu {NU};\n")
    _w(f"{d}/constant/turbulenceProperties", "dictionary", "turbulenceProperties",
       "simulationType laminar;\n")
    _w(f"{d}/constant/momentumTransport", "dictionary", "momentumTransport",
       "simulationType laminar;\n")
    OF._sh(["blockMesh"], d)
    inlet = _inlet_value(_inlet_profile(d))
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
    open(os.path.join(d, INLET_OK), "w").write("analytic, written by make_mesh\n")
    open(os.path.join(d, "MESHED"), "w").write("ok\n")
    return d


INLET_OK = "INLET_ANALYTIC"
VEC = re.compile(r"\(([-+\d.eE]+) [-+\d.eE]+ [-+\d.eE]+\)")


def _inlet_profile(d):
    """The analytic inlet velocity (DFG 2D-2, U_m = 1.5) at the inlet face centres."""
    _, cy = _patch_centres(d, "inlet")
    return 4.0 * UM2 * cy * (H - cy) / H**2


def _inlet_value(prof):
    return "nonuniform List<vector>\n%d\n(\n%s\n)" % (
        len(prof), "\n".join(f"({u:.15g} 0 0)" for u in prof))


def _inlet_block(s):
    """(start, end) of the body of the inlet entry in the boundaryField of a field file."""
    b0 = s.index("boundaryField")
    m = re.search(r"\binlet\b\s*\{", s[b0:])
    j = b0 + m.end()
    return j, s.index("}", j)


def _inlet_deviation(d, t="0"):
    """Largest |u_x - analytic| over the inlet faces of <d>/<t>/U (inf if not a list)."""
    s = open(os.path.join(d, t, "U")).read()
    j, k = _inlet_block(s)
    v = np.array([float(x) for x in VEC.findall(s[j:k])])
    prof = _inlet_profile(d)
    return float(np.max(np.abs(v - prof))) if len(v) == len(prof) else float("inf")


def _restore_inlet(d, t="0"):
    """mapFields -consistent maps the values of every patch, the fixedValue inlet
    included: write the analytic inlet of make_mesh back into <d>/<t>/U and check it."""
    f = os.path.join(d, t, "U")
    s = open(f).read()
    before = _inlet_deviation(d, t)
    j, k = _inlet_block(s)
    s = s[:j] + f"\n        type fixedValue;\n        value {_inlet_value(_inlet_profile(d))};\n    " + s[k:]
    open(f, "w").write(s)
    after = _inlet_deviation(d, t)
    if not after < 1e-12:
        raise RuntimeError(f"the analytic inlet could not be restored in {f} (deviation {after:.2e})")
    open(os.path.join(d, INLET_OK), "w").write(
        f"analytic, restored after mapFields (largest change {before:.3e})\n")
    return before


def _check_inlet(na):
    """Refuse a mapped mesh whose inlet was not restored (spun up before the fix)."""
    d = mesh_dir(na)
    if os.path.exists(os.path.join(d, "MAPPED")) and not os.path.exists(os.path.join(d, INLET_OK)):
        raise RuntimeError(f"mesh {na} was mapped before the inlet fix (no {INLET_OK}): its inlet is the "
                           f"coarser mesh's profile. Move {d} aside and spin it up again.")


def spun(na):
    f = os.path.join(mesh_dir(na), "SPUN_UP")
    return json.load(open(f)) if os.path.exists(f) else None


COURANT = re.compile(r"^Courant Number mean: ([-+\d.eE]+) max: ([-+\d.eE]+)", re.M)


def courant(logs):
    """Largest Courant number over the log files, and in the last 100 steps."""
    mx = []
    for f in logs:
        if os.path.exists(f):
            mx += [float(m.group(2)) for m in COURANT.finditer(open(f, errors="replace").read())]
    if not mx:
        return None, None
    return float(max(mx)), float(max(mx[-100:]))


def _stepped(log):
    """True if the solver log records at least one time step."""
    return os.path.exists(log) and re.search(r"^Time = ", open(log, errors="replace").read(),
                                                  re.M) is not None


def _clean_times(d):
    """Back to the initial state in 0/ (from rest, or mapped): remove later time
    directories and the decomposition of a failed attempt."""
    for x in os.listdir(d):
        if x.startswith("processor") or (re.fullmatch(r"[\d.eE+-]+", x) and x != "0"):
            shutil.rmtree(os.path.join(d, x), ignore_errors=True)


def _drop_newest_checkpoint(d):
    """Before resuming an interrupted spin-up: the newest checkpoint may have been cut
    off mid-write (a SLURM time-out), so it is removed, in the case and in every
    processor directory; the run resumes from the one before."""
    dirs = [d] + sorted(glob.glob(os.path.join(d, "processor*")))
    times = {x for q in dirs if os.path.isdir(q) for x in os.listdir(q)
             if re.fullmatch(r"[\d.eE+-]+", x) and float(x) > 0}
    if times:
        newest = max(times, key=float)
        for q in dirs:
            shutil.rmtree(os.path.join(q, newest), ignore_errors=True)


def spin(na, prev=None):
    """Spin up mesh na (from rest, or mapped from the spun-up mesh `prev`).

    The attempt number is kept in SPIN_ATTEMPT, so a spin-up resubmitted after a
    time-out resumes at the time step it was using; each attempt and each resumed
    segment has its own log (log.spin_a<attempt>_s<segment>)."""
    if spun(na):
        _check_inlet(na)
        return spun(na)
    d = make_mesh(na)
    n = NPROCS[na]                 # the smoke test too, so that it covers the parallel path
    if prev is None:
        t_end = T_SPIN0
    else:
        src = mesh_dir(prev)
        if not os.path.exists(os.path.join(d, "MAPPED")):
            if os.path.exists(os.path.join(d, INLET_OK)):
                os.remove(os.path.join(d, INLET_OK))
            OF._sh(["mapFields", src, "-consistent", "-sourceTime", "latestTime"], d)
            change = _restore_inlet(d)
            print(f"n_arc {na}: mapped from n_arc {prev}; analytic inlet restored "
                  f"(mapFields had changed it by up to {change:.3e})", flush=True)
            open(os.path.join(d, "MAPPED"), "w").write(str(prev))
        _check_inlet(na)
        t_end = T_SPIN_MAPPED
    fa = os.path.join(d, "SPIN_ATTEMPT")
    attempt = int(open(fa).read().strip() or 0) if os.path.exists(fa) else 0
    attempt = min(attempt, SPIN_RETRIES)
    err = None
    while attempt <= SPIN_RETRIES:
        dt = SPIN_DT[na] / 2 ** attempt
        open(fa, "w").write(f"{attempt}\n")
        _w(f"{d}/system/fvSolution", "dictionary", "fvSolution", fv_solution(1e-4, 1e-8, 1))
        _w(f"{d}/system/controlDict", "dictionary", "controlDict",
           OF.control_dict(0, t_end, dt, t_end / 5))
        seg = len(glob.glob(os.path.join(d, f"log.spin_a{attempt}_s*")))
        if seg:
            _drop_newest_checkpoint(d)                 # resuming after a time-out
        log = os.path.join(d, f"log.spin_a{attempt}_s{seg}")
        try:
            OF.run("pimpleFoam", d, nprocs=n, log=log)
            err = None
            break
        except RuntimeError as exc:
            if not _stepped(log):
                # decomposePar, mpirun or the solver failed before the first time step: an
                # installation or launch problem, not instability, so the time step is kept
                # (on Hamilton such failures used up both halvings of revision 4)
                raise RuntimeError(f"spin-up of n_arc {na} did not start (no time step; the "
                                   f"time step is not reduced): {str(exc)[-700:]}")
            err = exc
            print(f"spin-up of n_arc {na} at dt {dt:g} failed (Courant max "
                  f"{courant([log])[0]}); restarting from the initial state at dt {dt / 2:g}",
                  flush=True)
            _clean_times(d)
            attempt += 1
    if err is not None:
        raise RuntimeError(f"spin-up of n_arc {na} failed down to dt {dt:g}: {str(err)[-500:]}")
    if n > 1:
        OF._sh(["reconstructPar", "-latestTime"], d)
        for p in [x for x in os.listdir(d) if x.startswith("processor")]:
            shutil.rmtree(os.path.join(d, p))
    co, co_last = courant(sorted(glob.glob(os.path.join(d, f"log.spin_a{attempt}_s*"))))
    mp = os.path.join(d, "MAPPED")
    info = dict(t0=float(OF.latest_time(d)), dt=dt,
                mapped_from=int(open(mp).read().strip()) if os.path.exists(mp) else None,
                courant_max=co, courant_max_last100=co_last, attempts=attempt + 1,
                inlet="analytic", inlet_deviation_end=_inlet_deviation(d, OF.latest_time(d)))
    json.dump(info, open(os.path.join(d, "SPUN_UP"), "w"))
    return info


def tid(th):
    return "na%d_dt%g_o%g_p%g" % (th["n_arc"], th["dt"], th["outer_tol"], th["p_tol"])


def TASKS():
    T = []
    for na in NA:
        for dt in DT:
            for o in OUTER:
                for p in PTOL:
                    th = dict(n_arc=na, dt=dt, outer_tol=o, p_tol=p)
                    T.append(dict(id=tid(th), theta=th, nprocs=NPROCS[na], n_arc=na))
    for th in REFS:
        T.append(dict(id="ref_" + tid(th), theta=th, nprocs=NPROCS[th["n_arc"]], n_arc=th["n_arc"]))
    return T


def run(task, out):
    th = task["theta"]
    na, dt = th["n_arc"], th["dt"]
    info = spun(na)
    if info is None:
        raise RuntimeError(f"mesh {na} is not spun up (run: python dfg2d2_hpc.py spin)")
    _check_inlet(na)
    t0 = info["t0"]
    md = mesh_dir(na)
    case = os.path.join(out, "case")
    n = task["nprocs"]
    if not os.path.exists(os.path.join(case, "STARTED")):
        shutil.rmtree(case, ignore_errors=True)
        os.makedirs(case)
        for sub in ("constant", "system"):
            shutil.copytree(os.path.join(md, sub), os.path.join(case, sub))
        shutil.copytree(os.path.join(md, OF.latest_time(md)), os.path.join(case, OF.latest_time(md)))
        open(os.path.join(case, "STARTED"), "w").write("ok\n")
    _w(f"{case}/system/fvSolution", "dictionary", "fvSolution",
       fv_solution(th["outer_tol"], th["p_tol"], 50))
    fields = os.environ.get("ERRBUDGET_FO", "1") == "0"     # testing without function objects
    fo = "" if fields else OF.force_coeffs_fo("fc", magU=1.0, lref=0.1, aref=0.1 * 0.01, every=1)
    _w(f"{case}/system/controlDict", "dictionary", "controlDict",
       OF.control_dict(t0, t0 + T_TOT, dt, dt if fields else T_TOT / 8, functions=fo,
                       write_control="timeStep" if fields else "runTime")
       .replace(f"writeInterval {dt:.10g};", "writeInterval 1;") if fields else
       OF.control_dict(t0, t0 + T_TOT, dt, T_TOT / 8, functions=fo))
    seg = len([x for x in os.listdir(case) if x.startswith("log.run")])
    log = os.path.join(case, f"log.run{seg}")
    failed = None
    try:
        OF.run("pimpleFoam", case, nprocs=n, log=log)
    except RuntimeError as exc:
        failed = str(exc)[-800:]
    # cost and convergence from every log segment (a resumed segment repeats the
    # steps after the last checkpoint: the later record of each step is kept, and
    # CPU time is accumulated per segment)
    steps, offset = {}, 0.0
    for k in range(seg + 1):
        S = OF.log_steps(os.path.join(case, f"log.run{k}"))
        nxt = [s for s in S]
        for tt, ex, its, nc in nxt:
            steps[round(tt, 10)] = (ex + offset, its, nc)
        if S:
            offset += S[-1][1]
    ts = np.array(sorted(steps))
    cpu = np.array([steps[t][0] for t in ts]) * n          # core-seconds
    its = np.array([steps[t][1] for t in ts])
    nc = np.array([steps[t][2] for t in ts])
    co, co_last = courant([os.path.join(case, f"log.run{k}") for k in range(seg + 1)])
    res = dict(theta=th, nprocs=n, t0=t0, failed=failed, steps=len(ts),
               nonconverged_steps=int(nc.sum()) if len(nc) else 0,
               max_outer=int(its.max()) if len(its) else None,
               cpu_total=float(cpu[-1]) if len(cpu) else None,
               courant_max=co, courant_max_last100=co_last, inlet="analytic",
               spun_up=info)
    try:
        if fields:
            from openfoam_cylinder import _wall_geometry
            from openfoam_dfg2d2 import _coeffs
            if n > 1:
                OF._sh(["reconstructPar"], case)
            geo = _wall_geometry(case, "cylinder")
            tt = sorted((float(x), x) for x in os.listdir(case)
                        if re.fullmatch(r"[\d.eE+-]+", x) and float(x) > t0 + 1e-9)
            if not tt:
                raise RuntimeError("no fields written")
            C = np.array([_coeffs(os.path.join(case, x), geo) for _, x in tt])
            t, cd, cl = np.array([v for v, _ in tt]), C[:, 0], C[:, 1]
        else:
            t, cd, cl = OF.read_coefficients(case, "fc")
        np.savez_compressed(os.path.join(out, "series.npz"), t=t, cd=cd, cl=cl, ts=ts,
                            cpu=cpu, its=its, nc=nc)
        res["series_points"] = int(len(t))
        res["t_last"] = float(t[-1])
    except RuntimeError as exc:
        res["no_series"] = str(exc)
        if len(ts) == 0:
            # nothing ran: an installation problem, not a property of the setting;
            # raise so that the task is retried after the problem is fixed
            raise RuntimeError(f"pimpleFoam did not start:\n{failed}")
    if failed is None or "series_points" in res:
        # keep only what the analysis needs
        for x in os.listdir(case):
            if x.startswith("processor") or re.fullmatch(r"[\d.eE+-]+", x):
                shutil.rmtree(os.path.join(case, x), ignore_errors=True)
    json.dump(res, open(os.path.join(out, "result.json"), "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "spin":
        ALL = sorted(set(NA + [r["n_arc"] for r in REFS]))
        levels = [int(x) for x in sys.argv[2:]] or ALL
        failed = []
        for na in sorted(levels):
            # map from the next coarser mesh that is spun up with the analytic inlet, if any
            # (a level whose spin-up failed is skipped: the next one maps from the last that
            # succeeded); submit.sh runs one level per job, coarse to fine
            coarser = [m for m in ALL if m < na and spun(m) and
                       (not os.path.exists(os.path.join(mesh_dir(m), "MAPPED")) or
                        os.path.exists(os.path.join(mesh_dir(m), INLET_OK)))]
            prev = coarser[-1] if coarser else None
            mf = os.path.join(mesh_dir(na), "SPIN_FAILED")
            try:
                print("spin", na, "from", prev, spin(na, prev), flush=True)
                if os.path.exists(mf):
                    os.remove(mf)
            except Exception as exc:                      # noqa: BLE001 - report, go on
                failed.append(na)
                os.makedirs(mesh_dir(na), exist_ok=True)
                open(mf, "w").write(f"{exc}\n")
                print(f"SPIN-UP FAILED for n_arc {na}: {exc}", flush=True)
        # the summary file lists every level whose last spin-up failed
        pf = os.path.join(ROOT, "SPIN_FAILED")
        bad = [m for m in ALL if os.path.exists(os.path.join(mesh_dir(m), "SPIN_FAILED"))]
        if bad:
            open(pf, "w").write(" ".join(map(str, bad)) + "\n")
        elif os.path.exists(pf):
            os.remove(pf)
        if failed:
            print(f"SPIN-UP FAILED for n_arc {failed}; their runs will not start", flush=True)
        for na in sorted(levels):
            print("spun up:", na, spun(na), flush=True)
        sys.exit(1 if failed else 0)         # a failed level stops the jobs that depend on it
