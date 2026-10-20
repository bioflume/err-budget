"""J1: DFG 2D-1 with a current OpenFOAM release (referee points 7 and 26 of revision 3;
revision 4: function object fixed and added to every run, stopping-criterion study).

Three sets of tasks:

  fo_na<N>    N = 8 ... 192 at the tightest solver settings (residualControl 1e-8,
              GAMG 1e-10, 20,000 iterations at most): c_D and c_L from OpenFOAM's
              forceCoeffs function object and from the article's post-processing
              of the written fields (exact cell centroids, and the vertex-mean
              centres of version 1.1), so that the two can be compared and the
              convergence followed beyond n_arc = 64.  Each task also records the
              wall distance from both cell centres and the viscous drag computed
              with each (the size of the version-1.1 bias at every mesh level).
  grid_<idx>  the 840 settings of the article's grid with the current release,
              through the article's driver (openfoam_cylinder.DFGCylinder.run,
              unchanged), plus the forceCoeffs function object (obs cd_fo, cl_fo)
              and, from the solver log, why every linear solve stopped (info "stop").
  rt_<...>    the stopping-criterion study at n_arc = 32, nCellsInCoarsestLevel 10:
              every residualControl x GAMG tolerance of the grid with relTol 0 and
              0.01 for both linear solvers (relTol 0.1, the article's value, is the
              grid itself).  With relTol 0 only the absolute tolerances stop the
              linear solves, so the change of the answer from relTol 0.1 shows
              whether the absolute tolerance was active.

Revision 4 corrections.  The function object of the fo tasks was written with
writeControl writeTime and writeInterval 20,000, i.e. at the 20,000th field write,
which never happens: it now writes every iteration (writeControl timeStep,
writeInterval 1) and the last line is the converged value.

    python -m farm dfg2d1_hpc array K NK [nprocs=1|8]
    ERRBUDGET_SMOKE=1  a few small tasks of each kind (submit.sh smoke)
    ERRBUDGET_FO=0     no function objects (testing)
"""
from __future__ import annotations

import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
import ofmodern as OF
import openfoam_cylinder as OC
from openfoam_cylinder import (DFGCylinder, REF, UBAR, D, DEPTH, NU, fv_solution, _w,
                               _force_coeffs, _wall_geometry, _internal_vectors)
from exp46_dfg_build import LAD, NAMES

FO_LEVELS = [8, 12, 16, 24, 32, 48, 64, 96, 128, 192]
RT_NARC, RT_NC = 32, 10
RT_RC = list(LAD[1])                 # residualControl of the grid, 1e-2 ... 1e-7
RT_TOL = list(LAD[2])                # GAMG tolerance of the grid, 1e-4 ... 1e-8
RT_REL = [0.0, 0.01]                 # relTol 0.1 is the grid (n_arc 32, n_coarsest 10)
GRID = list(itertools.product(*[range(len(l)) for l in LAD]))
SMOKE = os.environ.get("ERRBUDGET_SMOKE") == "1"
if SMOKE:
    FO_LEVELS = [8, 16]
    RT_NARC, RT_RC, RT_TOL, RT_REL = 8, [1e-3], [1e-6], [0.0]
    GRID = [(0, 1, 2, 0), (1, 3, 4, 0)]
FO_ON = os.environ.get("ERRBUDGET_FO", "1") != "0"
_P = {}


def _fo_block():
    return OF.force_coeffs_fo("fc", magU=UBAR, lref=D, aref=D * DEPTH,
                              cofr=(0.2, 0.2, DEPTH / 2), every=1)


# --- the article's driver, used unchanged, with three additions -------------------
# (1) the forceCoeffs function object in the controlDict of every run;
# (2) relTol of both linear solvers settable (0.1 in the article);
# (3) the solver's standard output kept, for the stopping-criterion counts.
# openfoam_cylinder looks control_dict, fv_solution and subprocess up as module
# globals at call time, so replacing them here changes nothing else in it.
_CD0, _FV0 = OC.control_dict, OC.fv_solution
_REL = {"v": 0.1}


def _control_dict_fo(end=OC.MAX_ITER):
    s = _CD0(end)
    return s + "functions\n{\n" + _fo_block() + "}\n" if FO_ON else s


def _fv_solution_rel(*a, **k):
    s = _FV0(*a, **k)
    assert s.count("relTol 0.1;") == 2, "fv_solution of openfoam_cylinder changed"
    return s if _REL["v"] == 0.1 else s.replace("relTol 0.1;", f"relTol {_REL['v']:g};")


class _Tap:
    """subprocess for openfoam_cylinder that remembers the last standard output."""
    last = None

    def run(self, *a, **k):
        r = subprocess.run(*a, **k)
        _Tap.last = r.stdout
        return r

    def __getattr__(self, name):
        return getattr(subprocess, name)


OC.control_dict, OC.fv_solution, OC.subprocess = _control_dict_fo, _fv_solution_rel, _Tap()

SOLVE = re.compile(r"^\s*(\w+):\s+Solving for (\w+), Initial residual = ([-+\d.eE]+), "
                   r"Final residual = ([-+\d.eE]+), No Iterations (\d+)", re.M)


def parse_solves(log):
    """Per outer iteration, the linear solves of a simpleFoam log:
    list of [(field, initial, final, iterations), ...]."""
    its = []
    for block in re.split(r"^Time = ", log, flags=re.M)[1:]:
        its.append([(m.group(2), float(m.group(3)), float(m.group(4)), int(m.group(5)))
                    for m in SOLVE.finditer(block)])
    return its


def classify(init, final, n, tol, rel, max_iter=1000):
    """Why a linear solve stopped (OpenFOAM: converged when final < tolerance or,
    if relTol > 0, final < relTol * initial; no iteration if initial < tolerance).
    Residuals are printed to six digits, so a solve within that of both limits
    counts as 'both'."""
    if n == 0:
        return "zero"
    a = final < tol * (1 + 1e-5)
    r = rel > 0 and final < rel * init * (1 + 1e-5)
    if a and r:
        return "both"
    if a:
        return "abs"
    if r:
        return "rel"
    return "max" if n >= max_iter else "other"


def stop_counts(its, tol_p, rel, last=20):
    """Counts of the stopping reasons of the p and U solves over the whole run and
    over its last `last` outer iterations (the regime that sets the answer)."""
    def count(sel):
        out = {"p": {}, "U": {}}
        for it in sel:
            for f, i0, i1, n in it:
                g = "p" if f == "p" else ("U" if f.startswith("U") else None)
                if g is None:
                    continue
                c = classify(i0, i1, n, tol_p if g == "p" else 10 * tol_p, rel)
                out[g][c] = out[g].get(c, 0) + 1
        return out
    return dict(all=count(its), last=count(its[-last:]), last_n=min(last, len(its)),
                iters_logged=len(its))


def history(its):
    """Initial residuals and iteration counts per outer iteration (first solve of
    each field in the iteration)."""
    h = {}
    for k, it in enumerate(its):
        seen = set()
        for f, i0, i1, n in it:
            if f in seen:
                continue
            seen.add(f)
            h.setdefault(f, []).append((k + 1, i0, i1, n))
    return {f: np.array(v) for f, v in h.items()}


def wall_diagnostic(case, tdir=None):
    """Wall distance of the wall-cell centres (centroid and vertex mean) and, given
    the fields, the viscous drag computed with each."""
    geo = _wall_geometry(case, "cylinder")
    nh = geo["sf"] / np.linalg.norm(geo["sf"], axis=1)[:, None]
    dc = np.abs(np.einsum("ij,ij->i", geo["cc"] - geo["cf"], nh))
    dv = np.abs(np.einsum("ij,ij->i", geo["cc_vm"] - geo["cf"], nh))
    off = (dv - dc) / dc
    out = dict(n_faces=int(len(dc)), d_centroid_median=float(np.median(dc)),
               d_centroid_min=float(dc.min()), d_vertexmean_median=float(np.median(dv)),
               rel_offset_min=float(off.min()), rel_offset_median=float(np.median(off)),
               rel_offset_max=float(off.max()))
    if tdir is not None:
        U = _internal_vectors(os.path.join(tdir, "U"))[geo["cells"]]
        mag = np.linalg.norm(geo["sf"], axis=1)
        ut = U - np.einsum("ij,ij->i", U, nh)[:, None] * nh
        q = 0.5 * UBAR**2 * D * DEPTH
        for name, d in (("centroid", dc), ("vertexmean", dv)):
            F = (NU * mag[:, None] * ut / d[:, None]).sum(axis=0)
            out[f"cd_viscous_{name}"] = float(F[0] / q)
    return geo, out


def TASKS():
    T = [dict(id=f"fo_na{n}", kind="fo", n_arc=n, nprocs=8 if n >= 96 else 1)
         for n in FO_LEVELS]
    for idx in GRID:
        T.append(dict(id="grid_" + "_".join(map(str, idx)), kind="grid", idx=list(idx),
                      nprocs=1))
    for rel in RT_REL:
        for rc in RT_RC:
            for t in RT_TOL:
                th = dict(n_arc=RT_NARC, res_ctrl=rc, lin_tol=t, n_coarsest=RT_NC)
                T.append(dict(id=f"rt_na{RT_NARC}_rc{rc:g}_t{t:g}_rel{rel:g}", kind="rt",
                              theta=th, rel=rel, nprocs=1))
    return T


def _solver():
    if "p" not in _P:
        _P["p"] = DFGCylinder(workdir=tempfile.mkdtemp(prefix="j1_", dir=os.environ.get("TMPDIR")))
    return _P["p"]


def _run_setting(th, rel, out, keep_history=False):
    """One run through the article's driver; adds the function-object forces and the
    stopping-criterion counts."""
    p = _solver()
    _REL["v"] = rel
    _Tap.last = None
    try:
        obs, cost, wall, info = p.run(th)
    finally:
        _REL["v"] = 0.1
    info = dict(info)
    obs = dict(obs)
    log = _Tap.last or ""
    its = parse_solves(log)
    info["stop"] = stop_counts(its, float(th["lin_tol"]), rel)
    info["relTol"] = rel
    if FO_ON:
        try:
            t, cd, cl = OF.read_coefficients(os.path.join(p.root, "run"), "fc")
            obs["cd_fo"], obs["cl_fo"] = float(cd[-1]), float(cl[-1])
            info["fo_time"] = float(t[-1])
        except (RuntimeError, ValueError) as exc:
            obs["cd_fo"] = obs["cl_fo"] = float("nan")
            info["fo_error"] = str(exc)[-300:]
    if keep_history and its:
        np.savez_compressed(os.path.join(out, "residuals.npz"),
                            **{k: v for k, v in history(its).items()})
    return obs, cost, wall, info


def run(task, out):
    if task["kind"] == "grid":
        idx = tuple(task["idx"])
        th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
        obs, cost, wall, info = _run_setting(th, 0.1, out)
        json.dump(dict(idx=list(idx), theta=th, obs=obs, cost=cost, wall=wall, info=info),
                  open(os.path.join(out, "result.json"), "w"), indent=1)
        return
    if task["kind"] == "rt":
        th = dict(task["theta"])
        obs, cost, wall, info = _run_setting(th, float(task["rel"]), out, keep_history=True)
        json.dump(dict(theta=th, rel=task["rel"], obs=obs, cost=cost, wall=wall, info=info),
                  open(os.path.join(out, "result.json"), "w"), indent=1)
        return
    n = task["n_arc"]
    p = _solver()
    md = p.mesh(n)[0]
    case = os.path.join(out, "case")
    shutil.rmtree(case, ignore_errors=True)
    shutil.copytree(md, case)
    _, geo_diag = wall_diagnostic(case)            # mesh only, before the solver runs
    json.dump(geo_diag, open(os.path.join(out, "wall_geometry.json"), "w"), indent=1)
    end = 20000
    _w(f"{case}/system/fvSolution", "dictionary", "fvSolution", fv_solution(1e-8, 1e-10, 0.9))
    fo = _fo_block() if FO_ON else ""
    _w(f"{case}/system/controlDict", "dictionary", "controlDict",
       OF.control_dict(0, end, 1, end, functions=fo, app="simpleFoam",
                       write_control="timeStep").replace("adjustTimeStep no;", ""))
    log = os.path.join(case, "log")
    nproc = task["nprocs"]
    OF.run("simpleFoam", case, nprocs=nproc, log=log)
    if nproc > 1:
        OF._sh(["reconstructPar", "-latestTime"], case)
    txt = open(log).read()
    its = len(re.findall(r"^Time = ", txt, flags=re.M))
    conv = "SIMPLE solution converged" in txt
    last = OF.latest_time(case)
    tdir = os.path.join(case, last)
    geo, wd = wall_diagnostic(case, tdir)
    cd_c, cl_c = _force_coeffs(tdir, geo)
    cd_v, cl_v = _force_coeffs(tdir, geo, centres="cc_vm")
    res = dict(n_arc=n, iters=its, converged=conv,
               cd_fields=cd_c, cl_fields=cl_c, cd_fields_vertexmean=cd_v, cl_fields_vertexmean=cl_v,
               wall=wd)
    try:
        t, cd, cl = OF.read_coefficients(case, "fc")
        res.update(cd_fo=float(cd[-1]), cl_fo=float(cl[-1]), fo_time=float(t[-1]))
    except (RuntimeError, ValueError) as exc:
        res["fo_error"] = str(exc)
    for k in ("cd_fo", "cd_fields", "cd_fields_vertexmean"):
        if k in res:
            res["err_" + k] = abs(res[k] - REF["cd"]) / REF["cd"]
    for k in ("cl_fo", "cl_fields", "cl_fields_vertexmean"):
        if k in res:
            res["err_" + k] = abs(res[k] - REF["cl"]) / REF["cl"]
    json.dump(res, open(os.path.join(out, "result.json"), "w"), indent=1)
    shutil.rmtree(case, ignore_errors=True)
