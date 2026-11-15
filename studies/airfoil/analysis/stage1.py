#!/usr/bin/env python3
"""Stage-1 analysis: reference, admission, grading, and setting (a) (PROTOCOL.md 5-7).

    python3 analysis/stage1.py --runs DIR [--config DIR] [--out DIR] [--starts-every K]
                               [--write-settings]

  1. accounts for the 540 grid runs and the 5 reference runs;
  2. the reference (d) at the hardest design, its uncertainty u = |Q(d) - Q(c)|/|Q(d)|,
     the one-step-back siblings and the observed mesh order from n_surf 256, 362, 512;
     u_adm (5.2) and the targets admitted (>= 10 u_adm; 3 % if 1 % is not);
  3. grades errbudget 1.2.1 and S, SJ (eta = eps, eps/3, eps/10) and Ref. from every
     converged start, as Section 6.3-6.4 of the article (--starts-every K: every K-th
     start only, for tests);
  4. the user's case: errbudget through its Python API (`tune`, relative target, default
     start = middle of each list) on the stored runs, at 1 % and 0.3 %: setting, error
     against the reference, runs, search cost and K;
  5. with --write-settings, config/stage2_settings.json with (a)-(d).
Needs the errbudget engine (the article's src/, $ERRBUDGET_SRC, or engine/ of the tarball).
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, engine, fmt, md_table, reference_value, study, use_config, write_json  # noqa: E402
from accounting import account, report  # noqa: E402
import replay  # noqa: E402

LABELS = {"program": "errbudget", "seq-1": "S, eta=eps", "seq-3": "S, eta=eps/3", "seq-10": "S, eta=eps/10",
          "coord-1": "SJ, eta=eps", "coord-3": "SJ, eta=eps/3", "coord-10": "SJ, eta=eps/10", "ref": "Ref."}


def tensor_of(store, runs):
    T = {}
    for idx, r in zip(study.GRID, runs):
        j = store.get(r)
        if j is None:
            T[idx] = (math.nan, math.inf, {"converged": False, "missing": True})
            continue
        ok = j["status"] == "ok" and j["converged"]
        T[idx] = (float(j["ld"]) if ok else math.nan, float(j["cost"]["cells_iters"]),
                  {"converged": ok, "core_s": float(j["cost"]["core_s"]), "stop": j["stop_reason"]})
    return T


def reference(store, refruns, cautious_run):
    """PROTOCOL.md 5.2 as amended by A5: the reference value carries the far-field correction
    Delta_R; the siblings and the mesh order use the uncorrected values (all at R = 5000 but
    the radius sibling); u_adm adds 20 % of the correction."""
    d, d_mesh, d_dy1, d_rad, d6, ff = refruns
    q = {k: (store.value(r) if store.usable(r) else math.nan)
         for k, r in (("d", d), ("d_mesh", d_mesh), ("d_dy1", d_dy1), ("d_radius", d_rad), ("d6", d6))}
    qc = store.value(cautious_run) if store.usable(cautious_run) else math.nan
    qd = q["d"]
    out = dict(values=q, cautious=qc, converged={k: store.usable(r) for k, r in
                                                 (("d", d), ("d_mesh", d_mesh), ("d_dy1", d_dy1),
                                                  ("d_radius", d_rad), ("d6", d6), ("far_field_500", ff))})
    qref, ffinfo = reference_value(store, d, ff)
    out["reference_value"], out["far_field"] = qref, ffinfo
    out["u"] = abs(qref - qc) / abs(qref)
    out["sibling_changes"] = {k: abs(qd - q[k]) / abs(qd) for k in ("d_mesh", "d_dy1", "d_radius")}
    a, b = q["d"] - q["d_mesh"], q["d6"] - q["d"]
    p = math.log(abs(a) / abs(b)) / math.log(math.sqrt(2.0)) if a and b and math.isfinite(a + b) else math.nan
    out["mesh_order"] = p
    out["mesh_monotone"] = bool(a * b > 0) if math.isfinite(a + b) else None
    rich = math.nan
    if math.isfinite(p) and p > 0:
        rich = abs(b) / (math.sqrt(2.0) ** p - 1.0) / abs(qd)
    out["richardson_mesh"] = rich
    use_rich = math.isfinite(p) and p < 2 and math.isfinite(rich)
    u_ff = study.FAR_FIELD_UNCERTAINTY * abs(ffinfo["delta"]) / abs(qref)
    out["u_far_field"] = u_ff
    out["u_adm"] = (max(out["u"], rich) if use_rich else out["u"]) + u_ff
    out["u_adm_rule"] = (("max(u, Richardson of n_surf 362-512)" if use_rich else "u")
                         + " + 0.2 |Delta_R| / |L/D(d)| (A5)")
    return out


def user_case(errbudget, store, runs, T, target, qref):
    """errbudget's Python API on the stored runs, as a user would call it."""
    lad = study.ladders_now()
    params = [errbudget.Parameter(n, v, kind=k) for n, v, k in zip(study.NAMES, lad, study.KINDS)]
    by_setting = {tuple(study.setting_of_ladders(idx, lad).values()): (idx, r) for idx, r in zip(study.GRID, runs)}
    used = {}

    def solver(**kw):
        idx, r = by_setting[tuple(kw[n] for n in study.NAMES)]
        used[idx] = r
        v, cost, info = T[idx]
        if not np.isfinite(cost):
            raise errbudget.RunFailed(f"no run for {r.tag}")
        return {"ld": v if np.isfinite(v) else 0.0, "converged": bool(info["converged"])}, cost

    res = errbudget.tune(solver, params, target=target, observable="ld", start=study.default_start(),
                         relative=True)
    out = dict(target=target, certified=bool(res.certified), setting=res.settings,
               error_estimate=float(res.error_estimate), eps_abs=float(res.target), limiting=res.limiting,
               budget_shares=res.budget_shares, wasted=res.wasted, runs_used=len(used),
               search_cells_iters=float(sum(T[i][1] for i in used if np.isfinite(T[i][1]))),
               search_core_h=float(sum(T[i][2].get("core_s", 0.0) for i in used)) / 3600.0)
    if res.certified:
        idx = next(i for i in study.GRID if study.setting_of_ladders(i, lad) == res.settings)
        out["index"] = list(idx)
        out["ld"] = T[idx][0]
        out["cost_cells_iters"] = T[idx][1]
        out["core_h"] = T[idx][2].get("core_s", math.nan) / 3600.0
        out["K"] = out["search_cells_iters"] / T[idx][1]
        out["error_vs_reference"] = abs(T[idx][0] - qref) / abs(qref)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(os.environ.get("AIRFOIL_WORK", "."), "runs"))
    ap.add_argument("--config")
    ap.add_argument("--out")
    ap.add_argument("--starts-every", type=int, default=1)
    ap.add_argument("--write-settings", action="store_true")
    a = ap.parse_args(argv)
    cfg_dir = use_config(a.config)
    out_dir = a.out or os.path.join(os.path.dirname(os.path.abspath(a.runs)), "results")
    eb = engine()
    store = Store(a.runs)
    runs, refruns = study.sweep_stage1(), study.sweep_stage1ref()
    design = study.hardest_design()
    lad = study.ladders_now()
    print(f"hardest design m = {design[0]}, t = {design[1]}; engine {getattr(eb, '__version__', '?')} from "
          f"{os.path.dirname(eb.__file__)}")
    acc = [account(store, "stage1"), account(store, "stage1ref")]
    print("== 1. Accounting\n" + "\n".join(report(x) for x in acc))

    T = tensor_of(store, runs)
    cautious_idx = tuple(len(v) - 1 for v in lad)
    cautious_run = runs[study.GRID.index(cautious_idx)]
    ref = reference(store, refruns, cautious_run)
    qref = ref["reference_value"]                       # A5: L/D(d) + Delta_R
    print("\n== 2. Reference (d) and its uncertainty")
    print(md_table(["run", "L/D", "converged"],
                   [[k, fmt(v, 8), ref["converged"][k]] for k, v in ref["values"].items()]
                   + [["(c) cautious", fmt(ref["cautious"], 8), store.usable(cautious_run)],
                      ["far field: n_surf 91 at R = 5000 / 500", f"{fmt(ref['far_field']['ld_5000'], 8)} / "
                       f"{fmt(ref['far_field']['ld_500'], 8)}", ref["converged"]["far_field_500"]],
                      ["reference value L/D(d) + Delta_R (A5)", fmt(qref, 8),
                       f"Delta_R = {fmt(ref['far_field']['delta'], 4)}"]]))
    print(f"u = |Q(d) - Q(c)|/|Q(d)| = {ref['u']:.3g}; one-step-back changes: "
          + ", ".join(f"{k} {v:.3g}" for k, v in ref["sibling_changes"].items())
          + f"; mesh order from 256/362/512: {ref['mesh_order']:.3g} (monotone {ref['mesh_monotone']}); "
          f"Richardson error of the 362 mesh {ref['richardson_mesh']:.3g}")
    print(f"u_adm = {ref['u_adm']:.3g} ({ref['u_adm_rule']})")
    targets = []
    for tg in study.TARGETS:
        targets.append((tg, tg >= 10 * ref["u_adm"]))
    if not targets[0][1]:
        targets.append((study.FALLBACK_TARGET, study.FALLBACK_TARGET >= 10 * ref["u_adm"]))
    for tg, adm in targets:
        print(f"target {tg:.1%}: {'ADMITTED' if adm else 'not admitted (graded as information only)'}")

    print("\n== 3. Grading from every converged start (PROTOCOL.md 6)")
    starts = [i for i in study.GRID if T[i][2]["converged"]][::a.starts_every]
    grades = {}
    for tg, adm in targets:
        eps = tg * abs(qref)
        res = replay.grade(T, lad, study.hs_of(lad), study.KINDS, qref, eps, eb.run_protocol, dict(eb.V121),
                           starts=starts)
        rows = replay.summarise(res)
        grades[f"{tg}"] = dict(target=tg, admitted=adm, eps=eps, n_starts=res["n_starts"],
                               feasible=res["feasible"], good=res["good"], optimum=res["optimum"],
                               c_opt=res["c_opt"], rows=rows,
                               misses=[(list(s), list(q), e) for s, q, e in res["procs"]["program"]["misses"][:30]])
        print(f"\ntarget {tg:.1%} ({'admitted' if adm else 'NOT admitted'}): eps {eps:.4g}; {res['n_starts']} starts; "
              f"{res['feasible']} of {res['good']} converged settings meet it; cheapest {res['optimum']}")
        print(md_table(["procedure", "reached", "missed", "declined (reachable / not)", "cost ratio median (worst)",
                        "runs median", "K median", "distinct settings", "worst miss / eps"],
                       [[LABELS[p], r["reached"], r["missed"],
                         f"{r['declined']} ({r['declined_reachable']} / {r['declined_unreachable']})",
                         f"{fmt(r['penalty_median'], 3)} ({fmt(r['penalty_worst'], 3)})", fmt(r["runs_median"], 3),
                         fmt(r["K_median"], 3), r["distinct_settings"], fmt(r["worst_miss"], 3)]
                        for p, r in rows.items()]))
        if rows["program"]["reasons"]:
            print("declines by reason: " + ", ".join(f"{k} {v}" for k, v in rows["program"]["reasons"].items()))

    print("\n== 4. The user's case: default start, Python API, relative target")
    cases = {}
    for tg in study.TARGETS:
        c = user_case(eb, store, runs, T, tg, qref)
        cases[f"{tg}"] = c
        if c["certified"]:
            print(f"{tg:.1%}: setting {c['setting']}; L/D {c['ld']:.6g}, error vs reference {c['error_vs_reference']:.3g} "
                  f"({'within' if c['error_vs_reference'] <= tg else 'OUTSIDE'} the target); "
                  f"{c['runs_used']} runs, search {c['search_core_h']:.1f} core-h, K = {c['K']:.0f}; "
                  f"limiting: {c['limiting']}")
        else:
            print(f"{tg:.1%}: DECLINED (estimate {c['error_estimate']:.3g}); {c['runs_used']} runs, "
                  f"search {c['search_core_h']:.1f} core-h")

    settings = dict(b=study.tutorial_setting(), c=study.cautious_setting(), d=study.reference_setting())
    primary = cases[f"{study.TARGETS[0]}"]
    if primary["certified"]:
        settings["a"] = {k: (float(v) if k != "n_surf" else int(v)) for k, v in primary["setting"].items()}
    summary = dict(design=list(design), accounting=acc, reference=ref, targets=targets, grades=grades,
                   user_case=cases, stage2_settings=settings, starts_every=a.starts_every)
    write_json(os.path.join(out_dir, "stage1_summary.json"), summary)
    print(f"\nsummary written to {os.path.join(out_dir, 'stage1_summary.json')}")
    if a.write_settings:
        write_json(os.path.join(cfg_dir, "stage2_settings.json"),
                   dict(note="written by analysis/stage1.py (PROTOCOL.md 7)", settings=settings,
                        a_from=dict(target=study.TARGETS[0], start=study.default_start(),
                                    certified=primary["certified"], runs=primary["runs_used"],
                                    K=primary.get("K"))))
        print(f"stage-2 settings written to {cfg_dir}/stage2_settings.json: {sorted(settings)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
