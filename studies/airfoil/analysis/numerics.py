#!/usr/bin/env python3
"""Summary of the numerics sweep (study.sweep_numerics), on Hamilton or the laptop:

    python3 analysis/numerics.py            # $AIRFOIL_WORK/runs, or --runs DIR

For each run: variant, design, mesh, outcome, iterations, the final largest residual and
its equation, the iterations to bring the largest residual below 1e-4 and 1e-5, C_l and
how much it still moves over the last fifth of the iterations, the inviscid C_l of the
section (panel method), C_d, L/D, where the upper and lower boundary layers separate
(first wall face with reversed flow, x/c), and where the solution still changes most."""
from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, fmt, md_table, study  # noqa: E402
from panel import inviscid_cl  # noqa: E402


def first_below(hist, level):
    for it, r in hist:
        if r < level:
            return it
    return None


def cl_spread(j, frac=0.2):
    h = j.get("history") or {}
    it, cl = h.get("it") or [], h.get("cl") or []
    if not it:
        return math.nan
    cut = it[-1] - frac * (it[-1] - it[0])
    tail = [c for i, c in zip(it, cl) if i >= cut]
    m = sum(tail) / len(tail)
    return (max(tail) - min(tail)) / abs(m) if m else math.nan


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(os.environ.get("AIRFOIL_WORK", "."), "runs"))
    a = ap.parse_args(argv)
    store = Store(a.runs)
    inv = {}
    rows = []
    for i, r in enumerate(study.sweep_numerics()):
        key = (r.m, r.t)
        if key not in inv:
            inv[key] = inviscid_cl(r.m, r.t)
        case = f"m{r.m * 100:g} t{r.t * 100:g} n{r.n_surf}"
        j = store.get(r)
        if j is None:
            errs = store.errors.get(r.tag, [])
            rows.append([i, r.num, case, "no result" + (f" ({len(errs)} errors)" if errs else "")]
                        + [""] * 9 + [fmt(inv[key], 3), ""])
            continue
        res = j.get("residuals_final") or {}
        worst = max(res, key=res.get) if res else "-"
        hist = j.get("residual_history") or []
        f = j.get("fields") or {}
        wall = f.get("wall") or {}
        sep = "/".join(fmt(wall.get(s, {}).get("separation_x"), 2) for s in ("upper", "lower")) if wall else "-"
        share = (f.get("change") or {}).get("share") or {}
        where = max(share, key=share.get) + f" ({share[max(share, key=share.get)]:.0%})" if share else "-"
        rows.append([i, r.num, case, f"{j['status']}/{j['stop_reason']}", j["iterations"],
                     f"{fmt(res.get(worst, math.nan), 2)} ({worst})",
                     fmt(first_below(hist, 1e-4)), fmt(first_below(hist, 1e-5)),
                     fmt(j.get("cl"), 4), fmt(cl_spread(j), 2), fmt(j.get("cd"), 4), fmt(j.get("ld"), 4),
                     sep, fmt(inv[key], 3), where])
    print(md_table(["#", "variant", "case", "outcome", "iterations", "largest final residual",
                    "it. to 1e-4", "it. to 1e-5", "C_l", "C_l spread (last 20 %)", "C_d", "L/D",
                    "separation x/c (upper/lower)", "inviscid C_l", "still changing most in"], rows))


if __name__ == "__main__":
    main()
