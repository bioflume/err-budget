#!/usr/bin/env python3
"""Stage-2 analysis: the optimisation with each setting, its cost, and the payback
(PROTOCOL.md 8-11).

    python3 analysis/stage2.py --runs DIR [--config DIR] [--stage1 FILE] [--out DIR]

For each setting of config/stage2_settings.json ((a) errbudget, (b) tutorial,
(c) cautious, (d) reference):
  * L/D at the 25 designs; a full quadratic in (m, t), scaled to [-1, 1]^2, fitted by
    least squares; its maximum over the box (dense evaluation); the best grid design;
  * against (d): |dm|, |dt|, predicted optimal L/D against the reference's, the loss
    (reference surface at its own optimum minus at this setting's optimum, relative),
    the best grid design, the largest error over the designs;
  * cost: cells x iterations and core-hours over the 25 designs; for (a) also the search
    (from the Stage-1 summary) and the total; the payback n* = K'/(mu - 1) against each
    other setting, where K' is the search in units of a mean run of (a) and mu the ratio
    of mean costs per design;
  * the reference uncertainty u = |Q(d) - Q(c)|/|Q(d)| at every design.
A design whose run did not converge is left out of that setting's fit and reported.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, fmt, md_table, reference_value, study, use_config, write_json  # noqa: E402
from accounting import account, report  # noqa: E402

# the design box scaled to [-1, 1]^2 (A4: m in [0, 0.04], t in [0.08, 0.12])
M0, MS = (min(study.CAMBERS) + max(study.CAMBERS)) / 2, (max(study.CAMBERS) - min(study.CAMBERS)) / 2
T0, TS = (min(study.THICKNESSES) + max(study.THICKNESSES)) / 2, (max(study.THICKNESSES) - min(study.THICKNESSES)) / 2


def features(m, t):
    x, y = (np.asarray(m) - M0) / MS, (np.asarray(t) - T0) / TS
    return np.stack([np.ones_like(x), x, y, x * x, x * y, y * y], axis=-1)


def fit(designs, values):
    """Least-squares quadratic; returns coefficients, RMS residual, R^2."""
    m, t = np.array([d[0] for d in designs]), np.array([d[1] for d in designs])
    A, v = features(m, t), np.asarray(values, float)
    c, *_ = np.linalg.lstsq(A, v, rcond=None)
    r = v - A @ c
    ss = float(np.sum((v - v.mean()) ** 2))
    return c, float(np.sqrt(np.mean(r ** 2))), (1.0 - float(np.sum(r ** 2)) / ss) if ss > 0 else float("nan")


def surface(c, m, t):
    return features(m, t) @ c


def optimum(c, n=401):
    """Maximum of the fitted quadratic over the closed design box."""
    m = np.linspace(min(study.CAMBERS), max(study.CAMBERS), n)
    t = np.linspace(min(study.THICKNESSES), max(study.THICKNESSES), n)
    M, Tt = np.meshgrid(m, t, indexing="ij")
    F = surface(c, M, Tt)
    i, j = np.unravel_index(np.argmax(F), F.shape)
    return float(m[i]), float(t[j]), float(F[i, j])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(os.environ.get("AIRFOIL_WORK", "."), "runs"))
    ap.add_argument("--config")
    ap.add_argument("--stage1", help="stage1_summary.json (for the search cost of (a))")
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    use_config(a.config)
    out_dir = a.out or os.path.join(os.path.dirname(os.path.abspath(a.runs)), "results")
    store = Store(a.runs)
    settings = study.stage2_settings()
    acc = account(store, "stage2")
    print("== 1. Accounting\n" + report(acc))

    per = {}
    for key in sorted(settings):
        s = settings[key]
        vals, cost, coreh, missing = {}, [], [], []
        for d in study.DESIGNS:
            # A5: setting (d) is the settled reference with its far-field run; its value is
            # L/D(d) + Delta_R and its cost includes the far-field run
            r = study.reference_run(d, s) if key == "d" else study.run_at(d, s)
            j = store.get(r)
            if j is not None:
                extra = store.get(study.far_field_run(d)) if key == "d" else None
                cost.append(j["cost"]["cells_iters"] + (extra["cost"]["cells_iters"] if extra else 0.0))
                coreh.append((j["cost"]["core_s"] + (extra["cost"]["core_s"] if extra else 0.0)) / 3600.0)
            if key == "d":
                v, _ = reference_value(store, r, study.far_field_run(d))
                if math.isfinite(v):
                    vals[d] = v
            elif store.usable(r):
                vals[d] = store.value(r)
            else:
                missing.append([d[0], d[1], None if j is None else j["stop_reason"]])
        per[key] = dict(setting=s, values=vals, unusable=missing,
                        cells_iters=float(np.sum(cost)) if cost else math.nan,
                        core_h=float(np.sum(coreh)) if coreh else math.nan,
                        n_cost=len(cost))
        if len(vals) >= 10:
            c, rms, r2 = fit(list(vals), list(vals.values()))
            mo, to, fo = optimum(c)
            best = max(vals, key=vals.get)
            per[key].update(coef=c.tolist(), rms=rms, r2=r2, opt=dict(m=mo, t=to, ld=fo),
                            best_design=list(best), best_ld=vals[best])

    print("\n== 2. Response surfaces and optima")
    ref = per.get("d")
    rows = []
    for key, p in per.items():
        if "opt" not in p:
            rows.append([key, len(p["values"]), "-", "-", "-", "-", "-", "-", "-", "-"])
            continue
        o = p["opt"]
        row = [key, len(p["values"]), f"{o['m'] * 100:.2f} %", f"{o['t'] * 100:.2f} %", fmt(o["ld"], 6),
               f"({p['best_design'][0] * 100:g}, {p['best_design'][1] * 100:g})", fmt(p["rms"], 3), fmt(p["r2"], 4)]
        if ref and "opt" in ref and key != "d":
            rc = np.array(ref["coef"])
            f_ref_opt = float(surface(rc, ref["opt"]["m"], ref["opt"]["t"]))
            f_ref_here = float(surface(rc, o["m"], o["t"]))
            loss = (f_ref_opt - f_ref_here) / f_ref_opt
            common = [d for d in p["values"] if d in ref["values"]]
            errs = [abs(p["values"][d] - ref["values"][d]) / abs(ref["values"][d]) for d in common]
            p["against_reference"] = dict(dm=abs(o["m"] - ref["opt"]["m"]), dt=abs(o["t"] - ref["opt"]["t"]),
                                          ld_opt_error=(o["ld"] - ref["opt"]["ld"]) / ref["opt"]["ld"],
                                          loss=loss, same_best_design=p["best_design"] == ref["best_design"],
                                          max_error=max(errs) if errs else math.nan,
                                          median_error=float(np.median(errs)) if errs else math.nan)
            g = p["against_reference"]
            row += [f"{g['dm'] * 100:.2f} / {g['dt'] * 100:.2f}", f"{g['loss']:.2e}"]
        else:
            row += ["-", "-"]
        rows.append(row)
    print(md_table(["setting", "designs fitted", "m*", "t*", "L/D* (fit)", "best grid design", "fit RMS", "R^2",
                    "|dm| / |dt| vs (d) [%]", "loss vs (d)"], rows))
    for key, p in per.items():
        if p["unusable"]:
            print(f"setting {key}: designs without a converged run: {p['unusable']}")

    # reference uncertainty at every design
    if "c" in per and "d" in per:
        u = {d: abs(per["d"]["values"][d] - per["c"]["values"][d]) / abs(per["d"]["values"][d])
             for d in per["d"]["values"] if d in per["c"]["values"]}
        if u:
            worst = max(u, key=u.get)
            print(f"reference uncertainty u = |Q(d) - Q(c)|/|Q(d)|: median {np.median(list(u.values())):.3g}, "
                  f"largest {u[worst]:.3g} at {worst}")
        per["d"]["u_by_design"] = {f"{d[0]},{d[1]}": v for d, v in u.items()}

    print("\n== 3. Cost and payback")
    search = None
    if a.stage1 and os.path.exists(a.stage1):
        s1 = json.load(open(a.stage1))
        uc = s1["user_case"].get(f"{study.TARGETS[0]}")
        if uc and uc.get("certified"):
            search = dict(cells_iters=uc["search_cells_iters"], core_h=uc["search_core_h"], runs=uc["runs_used"],
                          K_hardest=uc["K"])
    rows = []
    mean = {k: p["cells_iters"] / max(p["n_cost"], 1) for k, p in per.items()}
    for key, p in per.items():
        tot = p["cells_iters"] + (search["cells_iters"] if (key == "a" and search) else 0.0)
        toth = p["core_h"] + (search["core_h"] if (key == "a" and search) else 0.0)
        rows.append([key, fmt(p["cells_iters"], 4), fmt(p["core_h"], 4),
                     fmt(search["core_h"], 4) if key == "a" and search else "-", fmt(toth, 4),
                     fmt(mean[key] / mean["a"], 3) if "a" in mean else "-"])
        p["total_cells_iters"], p["total_core_h"] = tot, toth
    print(md_table(["setting", "25 runs: cells x iterations", "25 runs: core-h", "search: core-h", "total core-h",
                    "cost per design / (a)"], rows))
    payback = {}
    if "a" in per and search:
        Kp = search["cells_iters"] / mean["a"]
        for key in per:
            if key == "a":
                continue
            mu = mean[key] / mean["a"]
            n_star = Kp / (mu - 1.0) if mu > 1.0 else math.inf
            payback[key] = dict(mu=mu, K_prime=Kp, n_star=n_star, within_25=n_star <= 25)
            if math.isfinite(n_star):
                print(f"against ({key}): mu = {mu:.3g}, K' = {Kp:.1f} runs of (a): the search pays for itself "
                      f"after {n_star:.1f} runs ({'within' if n_star <= 25 else 'beyond'} the 25 of this optimisation)")
            else:
                print(f"against ({key}): mu = {mu:.3g}: ({key}) costs less per run than (a), so the search never "
                      f"pays for itself against it")
    elif "a" not in per:
        print("no setting (a): errbudget declined from the default start (PROTOCOL.md 11)")

    summary = dict(accounting=acc, per_setting={k: {kk: vv for kk, vv in p.items() if kk != "values"}
                                                | {"values": {f"{d[0]},{d[1]}": v for d, v in p["values"].items()}}
                                                for k, p in per.items()},
                   search=search, payback=payback)
    write_json(os.path.join(out_dir, "stage2_summary.json"), summary)
    print(f"\nsummary written to {os.path.join(out_dir, 'stage2_summary.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
