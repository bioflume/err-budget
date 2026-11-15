#!/usr/bin/env python3
"""Stage-0 analysis (PROTOCOL.md 3.2, 3.3, 4.1, 4.2, 12; RUNBOOK.md, step 6).

    python3 analysis/stage0.py --runs DIR [--config DIR] [--out DIR] [--write-config]

Prints, from the 30 Stage-0 runs:
  1. the accounting of the runs;
  2. the measured cost model: core-seconds per cell and iteration, iterations to meet
     residualControl against mesh level and tolerance, the parallel speed-up;
  3. y+ achieved and the first-cell height measured on the wall;
  4. the predicted core-hours of Stages 1 and 2 with this model, against the guards;
  5. the mechanical decisions: the hardest design (4.2), the far-field report (A5; 3.2 withdrawn)
     and the numerics check (3.3).
With --write-config it writes config/stage1.json into the configuration directory.  If
the numerics check (3.3) triggers it prints STOP: the fallback is withdrawn (A3), and
the study stops for an amendment.  Needs NumPy only, so it
runs on a Hamilton login node.
"""
from __future__ import annotations

import argparse
import datetime as dt
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, fmt, md_table, study, use_config, write_json  # noqa: E402
from accounting import account, report  # noqa: E402

from af import costmodel  # noqa: E402

FF_LIMIT = 0.001             # eps/10: the far-field error at R = 5000 is reported against it (A5)
TIE = 0.10                   # PROTOCOL.md 4.2: within 10 % of the largest delta


def runs_by_group():
    s = study.sweep_stage0()
    return dict(A=s[0:3], A1=s[3:4], B=s[4:6], C=s[6:14], D=s[14:26])


def hardest(store, C):
    """PROTOCOL.md 4.2 from the corner runs at n_surf 128 and 181."""
    rows = []
    for d in study.CORNERS:
        r128 = next(r for r in C if (r.m, r.t) == d and r.n_surf == 128)
        r181 = next(r for r in C if (r.m, r.t) == d and r.n_surf == 181)
        ok = store.usable(r128) and store.usable(r181)
        q1, q2 = store.value(r128), store.value(r181)
        delta = abs(q2 - q1) / abs(q2) if ok else math.inf
        rows.append(dict(design=list(d), ld_128=q1, ld_181=q2, delta=delta, converged=ok))
    dmax = max(r["delta"] for r in rows)
    if math.isinf(dmax):
        tied = [r for r in rows if math.isinf(r["delta"])]
    else:
        tied = [r for r in rows if r["delta"] >= (1 - TIE) * dmax]
    pick = max(tied, key=lambda r: (r["design"][1], r["design"][0]))
    return pick["design"], rows


def far_field(store, D):
    """A5 (the check of PROTOCOL.md 3.2 withdrawn): at each corner, the settled runs D at
    R = 50, 500, 5000 give A in L/D = Q + A/R through R = 500 and 5000, and the far-field error
    of R = 5000, Delta_R = (L/D(5000) - L/D(500)) / 9, reported against eps/10.  The radius
    list is never changed (factor 1)."""
    rows = []
    for d in study.CORNERS:
        rr = {r.radius: r for r in D if (r.m, r.t) == d}
        vals = {R: (store.value(r) if store.usable(r) else math.nan) for R, r in sorted(rr.items())}
        q5, q05 = vals.get(5000.0, math.nan), vals.get(500.0, math.nan)
        ok = math.isfinite(q5) and math.isfinite(q05)
        delta = study.far_field_correction(q5, q05) if ok else math.nan
        A = (q5 - q05) / (1 / 5000.0 - 1 / 500.0) if ok else math.nan
        rows.append(dict(design=list(d), ld=vals, delta_R=delta, delta_rel=abs(delta) / abs(q5) if ok else math.nan,
                         A_rel=A / q5 if ok else math.nan,
                         above_eps10=bool(ok and abs(delta) / abs(q5) > FF_LIMIT)))
    return 1, rows


def numerics_check(store, runs):
    """PROTOCOL.md 3.3 (A3): triggered if any run diverged or > 2 of the 17 runs at 1e-6 are
    unconverged; the study then stops (the fallback is withdrawn)."""
    diverged = [r.tag for r in runs if (store.get(r) or {}).get("status") == "failed"]
    fine = [r for r in runs if r.res_ctrl == 1e-6]
    unconv = [r.tag for r in fine if store.get(r) is not None and not store.get(r).get("converged")]
    fallback = bool(diverged) or len(unconv) > 2
    return fallback, dict(diverged=diverged, unconverged_at_1e6=unconv, n_at_1e6=len(fine))


def cost_model(store, runs):
    """Core-seconds per cell-iteration, and iterations against level and tolerance."""
    serial = [r for r in runs if r.np == 1 and store.usable(r)]
    tci = [store.get(r)["cost"]["core_s"] / store.get(r)["cost"]["cells_iters"] for r in serial]
    out = dict(t_cell_iter=dict(median=float(np.median(tci)), min=float(np.min(tci)), max=float(np.max(tci)),
                                p90=float(np.percentile(tci, 90)), n=len(tci)))
    # iterations at res_ctrl 1e-5 against n_surf: log I = log I128 + b log(n/128)
    pts = [(r.n_surf, store.get(r)["iterations"]) for r in serial if r.res_ctrl == 1e-5]
    n = np.array([p[0] for p in pts], float)
    it = np.array([p[1] for p in pts], float)
    if len(set(n)) >= 2:
        b, a = np.polyfit(np.log(n / 128.0), np.log(it), 1)
    else:
        b, a = 2 * math.log(study.ITER_GROWTH) / math.log(2), math.log(study.ITERS_MID[1e-5])
    growth = math.sqrt(2.0) ** b
    i5 = math.exp(a)
    # one decade of residualControl: A' (256, 1e-6) against A (256, 1e-5), and the D runs (91, 1e-6)
    ratios = []
    g = runs_by_group()
    a256 = next(r for r in g["A"] if r.n_surf == 256)
    if store.usable(a256) and store.usable(g["A1"][0]):
        ratios.append(store.get(g["A1"][0])["iterations"] / store.get(a256)["iterations"])
    d_it = [store.get(r)["iterations"] for r in g["D"] if store.usable(r)]
    if d_it:
        ratios.append(float(np.median(d_it)) / (i5 * growth ** (math.log(91 / 128.0) / math.log(math.sqrt(2)))))
    ratio = float(np.median(ratios)) if ratios else study.ITERS_MID[1e-6] / study.ITERS_MID[1e-5]
    # Stage 0 measures 1e-5 and 1e-6 only.  Other tolerances: the larger of the straight
    # line in log10(res_ctrl) through the two and the assumed shape (study.ITERS_MID)
    # scaled to the measured 1e-5 point; the loose tolerances cost little either way.
    step = i5 * (ratio - 1.0)
    iters = {}
    for rc in (1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
        line = i5 + step * (-5 - math.log10(rc))
        shape = i5 * study.ITERS_MID[rc] / study.ITERS_MID[1e-5]
        iters[rc] = max(50.0, line, shape) if rc != 1e-5 else i5
    out.update(iters_at_128=iters, iter_growth_per_level=growth, iter_exponent_b=float(b),
               decade_ratio=ratio, fit_points=len(pts))
    # parallel speed-up (run B)
    b1, b4 = g["B"]
    if store.usable(b1) and store.usable(b4):
        w1, w4 = store.get(b1)["cost"]["wall_s"], store.get(b4)["cost"]["wall_s"]
        out["parallel"] = dict(wall_1=w1, wall_4=w4, speedup=w1 / w4, efficiency=w1 / w4 / 4,
                               core_s_ratio=store.get(b4)["cost"]["core_s"] / store.get(b1)["cost"]["core_s"],
                               ld_difference=abs(store.value(b4) - store.value(b1)) / abs(store.value(b1)))
    return out


def predict(model, hardest_design, factor=1):
    """Core-hours of Stages 1 and 2 with the measured model (median and 90th percentile);
    a sequenced run (A4) counts its first run too."""
    lad = [list(v) for v in study.LADDERS]
    ref = dict(study.REFERENCE)

    def ch(runs, t):
        runs = list(runs) + [s for s in (study.sequence_run(r) for r in runs) if s is not None]
        return sum(costmodel.cells(r.n_surf, r.dy1, r.radius) *
                   study.predicted_iters(r, model["iters_at_128"], model["iter_growth_per_level"]) * t
                   for r in runs) / 3600.0

    grid = [study.run_at(hardest_design, {n: lad[k][j] for k, (n, j) in enumerate(zip(study.NAMES, idx))})
            for idx in study.GRID]
    refruns = study.stage1ref_runs(hardest_design)
    mid = {n: lad[k][len(lad[k]) // 2] for k, n in enumerate(study.NAMES)}
    settings = dict(b=dict(mid, **study.TUTORIAL_TOLERANCES), c={n: lad[k][-1] for k, n in enumerate(study.NAMES)},
                    d=ref, a_low=mid, a_high={n: lad[k][-1] for k, n in enumerate(study.NAMES)})
    out = {}
    for label, t in (("median", model["t_cell_iter"]["median"]), ("p90", model["t_cell_iter"]["p90"])):
        s2 = {k: ch(study.stage2_reference_runs(s) if k == "d" else [study.run_at(d, s) for d in study.DESIGNS], t)
              for k, s in settings.items()}
        out[label] = dict(stage1=ch(grid, t), stage1ref=ch(refruns, t), stage2_b=s2["b"], stage2_c=s2["c"],
                          stage2_d=s2["d"], stage2_a_range=[s2["a_low"], s2["a_high"]],
                          stage2=s2["b"] + s2["c"] + s2["d"] + s2["a_high"])
        out[label]["total_1_2"] = out[label]["stage1"] + out[label]["stage1ref"] + out[label]["stage2"]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(os.environ.get("AIRFOIL_WORK", "."), "runs"))
    ap.add_argument("--config")
    ap.add_argument("--out", default=None)
    ap.add_argument("--write-config", action="store_true")
    a = ap.parse_args(argv)
    cfg_dir = use_config(a.config)
    out_dir = a.out or os.path.join(os.path.dirname(os.path.abspath(a.runs)), "results")
    store = Store(a.runs)
    g = runs_by_group()
    allruns = study.sweep_stage0()
    acc = account(store, "stage0")
    print("== 1. Accounting\n" + report(acc))

    print("\n== 2. Cost model")
    model = cost_model(store, allruns)
    t = model["t_cell_iter"]
    print(f"core-seconds per cell and iteration: median {t['median']:.3g} (range {t['min']:.3g}-{t['max']:.3g}, "
          f"90th percentile {t['p90']:.3g}; {t['n']} serial runs); assumed {study.T_CELL_ITER:.3g}")
    print(f"iterations at n_surf 128: " + ", ".join(f"res_ctrl {k:.0e}: {v:.0f}" for k, v in model["iters_at_128"].items())
          + f"  (x{model['iter_growth_per_level']:.2f} per level; assumed x{study.ITER_GROWTH})")
    rows = []
    for r in allruns:
        j = store.get(r)
        if j:
            rows.append([r.tag, r.np, j["cells"], j["iterations"], j["stop_reason"], fmt(j["ld"], 6),
                         fmt(j["cost"]["core_s"] / 60, 3), fmt(j.get("yplus", {}).get("mean"), 3),
                         fmt(j.get("yplus", {}).get("max"), 3),
                         fmt(j.get("first_cell", {}).get("median", math.nan) / r.dy1, 3)])
    print(md_table(["run", "np", "cells", "iterations", "stop", "L/D", "core-min", "y+ mean", "y+ max",
                    "first cell / dy1"], rows))
    if "parallel" in model:
        p = model["parallel"]
        print(f"parallel (n_surf 181): 1 core {p['wall_1']:.0f} s, 4 cores {p['wall_4']:.0f} s: speed-up "
              f"{p['speedup']:.2f}, efficiency {p['efficiency']:.0%}; L/D differs by {p['ld_difference']:.1e}")
    forces_ok = all(math.isfinite(store.value(r)) for r in allruns if store.get(r) and store.get(r)["status"] == "ok")
    decay_ok = all(store.get(r)["openfoam"].get("decay_control_echoed") for r in allruns if store.get(r))
    print(f"forces written for every run: {forces_ok}; decayControl echoed by every run: {decay_ok}")

    print("\n== 3. Decisions (mechanical)")
    design, hrows = hardest(store, g["C"])
    print(md_table(["corner (m, t)", "L/D at 128", "L/D at 181", "relative change", "converged"],
                   [[tuple(r["design"]), fmt(r["ld_128"], 6), fmt(r["ld_181"], 6), fmt(r["delta"], 3), r["converged"]]
                    for r in hrows]))
    print(f"hardest design (PROTOCOL.md 4.2): m = {design[0]}, t = {design[1]}")
    factor, frows = far_field(store, g["D"])
    print(md_table(["corner", "L/D at R = 50, 500, 5000 (settled)", "Delta_R / (L/D) at 5000",
                    "A / (L/D) in L/D = Q + A/R", "above eps/10"],
                   [[tuple(r["design"]), ", ".join(fmt(v, 6) for v in r["ld"].values()),
                     fmt(r["delta_rel"], 3), fmt(r["A_rel"], 3), r["above_eps10"]] for r in frows]))
    print("far field (A5; the check of PROTOCOL.md 3.2 is withdrawn): radius list 50, 500, 5000 and "
          "reference radius 5000 unchanged; the reference value carries Delta_R")
    fallback, ninfo = numerics_check(store, allruns)
    print(f"numerics check (PROTOCOL.md 3.3): {'TRIGGERED' if fallback else 'passed'} "
          f"({len(ninfo['diverged'])} diverged; {len(ninfo['unconverged_at_1e6'])} of {ninfo['n_at_1e6']} "
          f"runs at 1e-6 unconverged)")

    print("\n== 4. Predicted core-hours with the measured model")
    pred = predict(model, tuple(design), factor)
    for label, p in pred.items():
        print(f"{label:>6} t: Stage 1 {p['stage1']:.0f} + reference {p['stage1ref']:.0f}; Stage 2 (b) "
              f"{p['stage2_b']:.0f}, (c) {p['stage2_c']:.0f}, (d) {p['stage2_d']:.0f}, (a) "
              f"{p['stage2_a_range'][0]:.0f}-{p['stage2_a_range'][1]:.0f}; Stages 1-2 at most {p['total_1_2']:.0f}")
    if fallback:
        print("STOP (PROTOCOL.md 3.3, A3): the numerics check triggered; the fallback is withdrawn. "
              "Report the runs above; Stage 1 waits for an amendment")
    stop = pred["median"]["total_1_2"] > 900 or fallback
    if pred["median"]["total_1_2"] > 900:
        print("STOP (PROTOCOL.md 12): Stages 1-2 predicted above 900 core-hours; revisit before Stage 1")
    s1 = pred["median"]["stage1"] + pred["median"]["stage1ref"]
    if s1 > 400 or pred["median"]["stage2"] > 500:
        print("WARNING: a stage is predicted above its guard (400 / 500 core-hours)")

    summary = dict(date=dt.date.today().isoformat(), accounting={k: v for k, v in acc.items()},
                   cost_model=model, hardest=design, corners=hrows, radius_factor=factor, far_field=frows,
                   numerics_check_triggered=fallback, numerics_check=ninfo, predicted=pred,
                   forces_written=forces_ok, decay_control_echoed=decay_ok, stop=stop)
    write_json(os.path.join(out_dir, "stage0_summary.json"), summary)
    print(f"\nsummary written to {os.path.join(out_dir, 'stage0_summary.json')}")
    if a.write_config:
        if not acc["complete"]:
            print("NOT writing config: Stage 0 is not complete (resubmit, then run this again)")
            return 1
        if fallback:
            print("NOT writing config: the numerics check triggered (PROTOCOL.md 3.3, A3)")
            return 1
        write_json(os.path.join(cfg_dir, "stage1.json"),
                   dict(note="written by analysis/stage0.py (PROTOCOL.md 3.2, 4.2)", hardest=design,
                        radius_factor=factor, corners=hrows, far_field=frows, date=summary["date"]))
        print(f"config written to {cfg_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
