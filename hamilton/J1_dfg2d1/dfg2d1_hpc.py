"""J1: DFG 2D-1 with a current OpenFOAM release (referee points 7 and 26).

Two sets of tasks:

  fo_na<N>    N = 8 ... 192 at the tightest solver settings (residualControl 1e-8,
              GAMG 1e-10, 20,000 iterations at most): c_D and c_L from OpenFOAM's
              forceCoeffs function object and from the article's post-processing
              of the written fields (exact cell centroids, and the vertex-mean
              centres of version 1.1), so that the two can be compared and the
              convergence followed beyond n_arc = 64;
  grid_<idx>  the 840 settings of the article's grid with the current release
              (post-processing of the written fields, as in the article), so the
              program's results can be compared between OpenFOAM v1912 and the
              current release.
    python -m farm dfg2d1_hpc array K NK
"""
from __future__ import annotations

import itertools
import json
import math
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
from openfoam_cylinder import (DFGCylinder, REF, UBAR, D, DEPTH, fv_solution, _w,
                               _force_coeffs, _wall_geometry, _internal_scalars, _interp_angle)
from exp46_dfg_build import LAD, NAMES

FO_LEVELS = [8, 12, 16, 24, 32, 48, 64, 96, 128, 192]
SMOKE = os.environ.get("KNOBS_SMOKE") == "1"
if SMOKE:
    FO_LEVELS = [8, 16]
_P = {}


def TASKS():
    T = [dict(id=f"fo_na{n}", kind="fo", n_arc=n, nprocs=8 if n >= 96 else 1)
         for n in FO_LEVELS]
    if not SMOKE:
        for idx in itertools.product(*[range(len(l)) for l in LAD]):
            T.append(dict(id="grid_" + "_".join(map(str, idx)), kind="grid", idx=list(idx),
                          nprocs=1))
    return T


def _solver():
    if "p" not in _P:
        _P["p"] = DFGCylinder(workdir=tempfile.mkdtemp(prefix="j1_", dir=os.environ.get("TMPDIR")))
    return _P["p"]


def run(task, out):
    if task["kind"] == "grid":
        idx = tuple(task["idx"])
        th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
        obs, cost, wall, info = _solver().run(th)
        json.dump(dict(idx=list(idx), theta=th, obs=obs, cost=cost, wall=wall, info=info),
                  open(os.path.join(out, "result.json"), "w"), indent=1)
        return
    n = task["n_arc"]
    p = _solver()
    md = p.mesh(n)[0]
    case = os.path.join(out, "case")
    shutil.rmtree(case, ignore_errors=True)
    shutil.copytree(md, case)
    end = 20000
    _w(f"{case}/system/fvSolution", "dictionary", "fvSolution", fv_solution(1e-8, 1e-10, 0.9))
    fo = OF.force_coeffs_fo("fc", magU=UBAR, lref=D, aref=D * DEPTH, cofr=(0.2, 0.2, DEPTH / 2),
                            every=end)
    fo = fo.replace("writeControl    timeStep;", "writeControl    writeTime;")
    if os.environ.get("KNOBS_FO", "1") == "0":       # testing without function objects
        fo = ""
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
    geo = _wall_geometry(case, "cylinder")
    cd_c, cl_c = _force_coeffs(tdir, geo)
    cd_v, cl_v = _force_coeffs(tdir, geo, centres="cc_vm")
    res = dict(n_arc=n, iters=its, converged=conv,
               cd_fields=cd_c, cl_fields=cl_c, cd_fields_vertexmean=cd_v, cl_fields_vertexmean=cl_v)
    try:
        t, cd, cl = OF.read_coefficients(case, "fc")
        res.update(cd_fo=float(cd[-1]), cl_fo=float(cl[-1]), fo_time=float(t[-1]))
    except RuntimeError as exc:
        res["fo_error"] = str(exc)
    for k in ("cd_fo", "cd_fields", "cd_fields_vertexmean"):
        if k in res:
            res["err_" + k] = abs(res[k] - REF["cd"]) / REF["cd"]
    json.dump(res, open(os.path.join(out, "result.json"), "w"), indent=1)
    shutil.rmtree(case, ignore_errors=True)
