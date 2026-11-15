#!/usr/bin/env python3
"""Summary of the diagnostic sweep (study.sweep_diag), printed on Hamilton:

    python3 analysis/diag.py

For each run: numerics, outcome, iterations, L/D, C_l, C_d, how much L/D still moves
over the last tenth of the iterations, the final residuals and y+."""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, fmt, md_table, study  # noqa: E402


def main():
    store = Store(os.path.join(os.environ.get("AIRFOIL_WORK", "."), "runs"))
    rows = []
    for i, r in enumerate(study.sweep_diag()):
        j = store.get(r)
        if j is None:
            errs = store.errors.get(r.tag, [])
            rows.append([i, f"m{r.m*100:g} t{r.t*100:g} n{r.n_surf} rc{r.res_ctrl:.0e}", r.num,
                         "no result" + (f" ({len(errs)} errors)" if errs else ""), "", "", "", "", "", "", ""])
            continue
        w = j.get("ld_window", {})
        spread = (w.get("max", float("nan")) - w.get("min", float("nan"))) / abs(j["ld"]) if j.get("ld") else float("nan")
        res = j.get("residuals_final", {})
        rows.append([i, f"m{r.m*100:g} t{r.t*100:g} n{r.n_surf} rc{r.res_ctrl:.0e}", r.num,
                     f"{j['status']}/{j['stop_reason']}", j["iterations"], fmt(j["ld"], 5), fmt(j["cl"], 4),
                     fmt(j["cd"], 4), fmt(spread, 2), fmt(max(res.values()) if res else float("nan"), 2),
                     fmt(j.get("yplus", {}).get("mean"), 3)])
    print(md_table(["#", "case", "numerics", "outcome", "iterations", "L/D", "C_l", "C_d",
                    "L/D spread (last 10 %)", "largest final residual", "y+ mean"], rows))


if __name__ == "__main__":
    main()
