#!/usr/bin/env python3
"""Synthetic run JSONs in the driver's format, to test the analyses without OpenFOAM.

    python3 analysis/synth.py --out DIR [--sweeps stage0 stage1 stage1ref stage2] [--scale 1.0]

L/D = L/D*(m, t) (1 + e_mesh + e_wall + e_far + e_res + e_lin), with
  L/D*(m, t) = 95 - 6000 (m - 0.04)^2 - 2500 (t - 0.10)^2           (optimum m 4 %, t 10 %)
  e_mesh = -scale a(m, t) (64 / n_surf)^2, a = 0.03 + 0.5 m            (second order; largest
           for the most cambered corners, so (4, 8) and (4, 12) are the hardest)
  e_wall = -scale 0.006 (dy1 / 4e-5)^2,  e_far = -scale far (50 / R),  e_res = 20 res_ctrl s,
  e_lin = 2 gamg_tol s  (s = +-1, fixed per setting: the iterative error has no sign)
Cost: the generator's cell count times the iteration model of study.py; core-seconds at
T_TRUE per cell-iteration.  Runs at res_ctrl 1e-6 with gamg_tol 1e-6 on n_surf 256 stall
(unconverged), and the coarsest mesh at R 50 with res_ctrl 1e-6 diverges, so the
analyses meet both kinds of failure.
"""
from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import study, write_json  # noqa: E402

from af import costmodel  # noqa: E402
from af.cgrid import MESH_VERSION  # noqa: E402

T_TRUE = 3.5e-6


def ld_true(m, t):
    return 95.0 - 6000.0 * (m - 0.04) ** 2 - 2500.0 * (t - 0.10) ** 2


def model(run, scale=1.0, far=0.025):
    sgn = 1.0 if int(hashlib.md5(run.tag.encode()).hexdigest(), 16) % 2 else -1.0
    a = 0.03 + 0.5 * run.m
    e = (-scale * a * (64.0 / run.n_surf) ** 2 - scale * 0.006 * (run.dy1 / 4e-5) ** 2
         - scale * far * (50.0 / run.radius) + 20.0 * run.res_ctrl * sgn + 2.0 * run.gamg_tol * sgn)
    return ld_true(run.m, run.t) * (1.0 + e)


def record(run, scale=1.0, far=0.025, stall=True, diverge=True):
    cells = costmodel.cells(run.n_surf, run.dy1, run.radius)
    iters = int(round(study.predicted_iters(run)))
    status, conv, stop = "ok", True, "converged"
    if stall and run.res_ctrl == 1e-6 and run.gamg_tol == 1e-6 and run.n_surf == 256:
        conv, stop, iters = False, "stalled", iters + study.STALL_ITER
    if diverge and run.n_surf == 64 and run.radius == 50.0 and run.res_ctrl == 1e-6:
        status, conv, stop, iters = "failed", False, "diverged", 400
    ld = model(run, scale, far) if status == "ok" else math.nan
    if run.settle and conv:
        stop = "settled"                                 # A5
    cd = 0.009 * (1 + 0.5 * (64.0 / run.n_surf) ** 2)
    cpu = cells * iters * T_TRUE
    wall = cpu / (3.4 if run.np == 4 else run.np)
    yp = 0.6 * run.dy1 / 1e-5
    return dict(schema=1, mesh_version=MESH_VERSION, setup_version=study.SETUP_VERSION, tag=run.tag, status=status, converged=conv, stop_reason=stop, run=run.to_dict(),
                numerics=study.numerics(), setting=run.setting(), np=run.np, iterations=iters,
                max_iter=run.max_iter, cells=cells, ld=ld, cd=cd, cl=ld * cd if status == "ok" else math.nan,
                cm=-0.1, ld_window=dict(n=100, min=ld, max=ld),
                mesh=dict(generator=dict(cells=cells), checkMesh=dict(cells=cells, max_non_orth=20.0)),
                yplus=dict(n_faces=2 * run.n_surf, min=0.1 * yp, mean=yp, p95=2 * yp, max=3 * yp, frac_gt_1=0.2),
                first_cell=dict(min=0.97 * run.dy1, median=run.dy1, max=1.01 * run.dy1, target=run.dy1),
                cost=dict(cells_iters=float(cells * iters), cpu_s=cpu, core_s=cpu * run.np, wall_s=wall),
                openfoam=dict(version="v2606", decay_control_echoed=True), segments=1, errors=[],
                sweeps=[], synthetic=True, **_sequence(run, scale, far))


def _sequence(run, scale, far):
    """A4/A5: the first run of a sequenced run, settled, with its L/D (the far-field correction
    of the reference uses it)."""
    seq = study.sequence_run(run)
    if seq is None:
        return {}
    cells = costmodel.cells(seq.n_surf, seq.dy1, seq.radius)
    iters = int(round(study.predicted_iters(seq)))
    return dict(sequence=dict(tag=seq.tag, iterations=iters, converged=True, stop="settled",
                              ld=model(seq, scale, far), cells=cells, cells_iters=float(cells * iters),
                              cpu_s=cells * iters * T_TRUE, wall_s=cells * iters * T_TRUE, residuals_final={}))


def write(out, sweeps, scale=1.0, far=0.025, **kw):
    os.makedirs(out, exist_ok=True)
    n = 0
    for s in sweeps:
        for run in study.sweep(s):
            p = os.path.join(out, f"{run.tag}.json")
            if not os.path.exists(p):
                write_json(p, record(run, scale, far, **kw))
                n += 1
    return n


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--sweeps", nargs="+", default=["stage0"])
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--far", type=float, default=0.025)
    a = ap.parse_args()
    print(write(a.out, a.sweeps, a.scale, a.far), "runs written")
