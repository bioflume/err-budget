#!/usr/bin/env python3
"""The cost estimate of COST.md, before Stage 0 (PROTOCOL.md 12).

    python3 analysis/cost_estimate.py [--t 3e-6] [--hardest 0.06 0.16]

Cells are exact (the mesh generator); the time per cell and iteration and the
iterations to convergence are assumptions (af/study.py: T_CELL_ITER, ITERS_MID,
ITER_GROWTH) until Stage 0 measures them.  Setting (a) is unknown before Stage 1: it
is bounded by the default start below and the cautious setting above.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import md_table, study  # noqa: E402

from af import costmodel  # noqa: E402


def core_h(runs, t):
    """Predicted core-hours; a sequenced run (A4) counts its first run too."""
    runs = list(runs) + [s for s in (study.sequence_run(r) for r in runs) if s is not None]
    return sum(costmodel.predicted_core_seconds(r, t) for r in runs) / 3600.0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--t", type=float, default=study.T_CELL_ITER)
    ap.add_argument("--hardest", type=float, nargs=2, default=[0.06, 0.16])
    a = ap.parse_args(argv)
    d = tuple(a.hardest)
    lad = study.LADDERS
    grid = [study.run_at(d, study.setting_of_ladders(i, lad)) for i in study.GRID]
    ref = dict(study.REFERENCE)
    refruns = study.stage1ref_runs(d)                   # A5: settled, with the far-field run
    mid = study.MID
    cautious = {n: v[-1] for n, v in zip(study.NAMES, lad)}
    tut = dict(mid, **study.TUTORIAL_TOLERANCES)
    rows, tot = [], 0.0
    items = [("Stage 0 (check job: 2 tiny runs)", study.sweep_check(), "test", "15 min"),
             ("Stage 0", study.sweep_stage0(), "shared", "3-12 h limits")]
    for n in lad[0]:
        sub = [r for r in grid if r.n_surf == n]
        items.append((f"Stage 1, n_surf {n} (108 runs)", sub, "shared", None))
    items.append((f"Stage 1, reference runs ({len(refruns)})", refruns, "shared", "36 h / 72 h limits"))
    for key, s in (("(a) errbudget, if = default start", mid), ("(a) errbudget, if = cautious", cautious),
                   ("(b) tutorial", tut), ("(c) cautious", cautious), ("(d) reference", ref)):
        rr = study.stage2_reference_runs(s) if key == "(d) reference" else [study.run_at(dd, s) for dd in study.DESIGNS]
        items.append((f"Stage 2, {key}", rr, "shared", None))
    for label, runs, part, note in items:
        h = core_h(runs, a.t)
        longest = max(costmodel.predicted_core_seconds(r, a.t) for r in runs) / 3600.0
        cells = [costmodel.cells(r.n_surf, r.dy1, r.radius) for r in runs]
        rows.append([label, len(runs), f"{min(cells):,}-{max(cells):,}", f"{h:.1f}", f"{2 * h:.1f}",
                     f"{longest:.2f}"])
    print(md_table(["item", "runs", "cells", "core-h (central)", "core-h (x2)", "longest run (h)"], rows))
    s0 = core_h(study.sweep_check() + study.sweep_stage0(), a.t)
    s1 = core_h(grid + refruns, a.t)
    s2 = core_h([study.run_at(dd, s) for s in (cautious, tut, cautious) for dd in study.DESIGNS]
                + study.stage2_reference_runs(ref), a.t)
    print(f"\nStage 0 {s0:.1f}, Stage 1 {s1:.1f}, Stage 2 {s2:.1f} (with (a) = cautious, the upper bound); "
          f"total {s0 + s1 + s2:.0f} core-hours central, {2 * (s0 + s1 + s2):.0f} with the x2 margin")
    return 0


if __name__ == "__main__":
    sys.exit(main())
