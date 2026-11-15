#!/usr/bin/env python3
"""Summary of a convergence map, on Hamilton or the laptop: the convergence map of A3
(study.sweep_convmap) or the confirmation map of A4 (study.sweep_confirm):

    python3 analysis/convmap.py --sweep confirm      # $AIRFOIL_WORK/runs, or --runs DIR

For each design: how many of the 45 Stage-1 meshes converged at residualControl 1e-6, and
a table of the meshes with the outcome, iterations, the largest final residual and L/D;
then the reference runs."""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, fmt, md_table, study  # noqa: E402


def outcome(j):
    if j is None:
        return "no result"
    if j.get("converged"):
        return "converged"
    return f"{j['status']}/{j['stop_reason']}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(os.environ.get("AIRFOIL_WORK", "."), "runs"))
    ap.add_argument("--sweep", default="convmap", choices=["convmap", "confirm"])
    a = ap.parse_args(argv)
    store = Store(a.runs)
    runs = study.sweep(a.sweep)
    designs = []
    for r in runs:
        if r.n_surf <= max(study.LADDERS[0]) and (r.m, r.t) not in designs:
            designs.append((r.m, r.t))
    for d in designs:
        mine = [r for r in runs if (r.m, r.t) == d and r.n_surf <= max(study.LADDERS[0])]
        res = [store.get(r) for r in mine]
        n_conv = sum(bool(j and j.get("converged")) for j in res)
        n_done = sum(j is not None for j in res)
        print(f"\n== m {d[0] * 100:g} %, t {d[1] * 100:g} %: converged {n_conv} of {len(mine)} "
              f"({n_done} with a result)")
        rows = []
        for r, j in zip(mine, res):
            fr = (j or {}).get("residuals_final") or {}
            worst = max(fr, key=fr.get) if fr else None
            rows.append([r.n_surf, f"{r.dy1:.0e}", f"{r.radius:g}", outcome(j),
                         (j or {}).get("iterations", ""),
                         f"{fmt(fr[worst], 2)} ({worst})" if worst else "",
                         fmt((j or {}).get("ld"), 5), fmt((j or {}).get("cl"), 4)])
        print(md_table(["n_surf", "dy1", "R", "outcome", "iterations", "largest final residual",
                        "L/D", "C_l"], rows))
    print("\n== reference setting")
    rows = []
    for r in runs:
        if r.n_surf > max(study.LADDERS[0]):
            j = store.get(r)
            fr = (j or {}).get("residuals_final") or {}
            seq = (j or {}).get("sequence") or {}
            first = f"{seq['stop']}, {seq['iterations']}" if seq else ""
            rows.append([f"m {r.m * 100:g} %, t {r.t * 100:g} %", outcome(j), (j or {}).get("iterations", ""),
                         fmt(max(fr.values()), 2) if fr else "", fmt((j or {}).get("ld"), 6), first])
    print(md_table(["design", "outcome", "iterations", "largest final residual", "L/D",
                    "first run (A4): stop, iterations"], rows))


if __name__ == "__main__":
    main()
