"""errbudget's side of the comparison, recomputed by the code that grades PRE.

1. Table S-3.1 (paper/tab_pooled.tex at 79edc9a, the 2D-1 rows): the pooled replay from all
   780 usable starts (pooled.replay), compared with the table and with data/pooled_v1.2.1.json.
2. Table 7 / exp72: reached, cost ratio, runs and K from the 608 converged starts, from
   results/arm1_rows.csv (the same replay that gave PRE its data), against data/baselines_v6.json.
3. errbudget-3p: the program on the 3-parameter problem at nCellsInCoarsestLevel = 10, from the
   164 converged starts there.

    python check_replay.py        # writes results/errbudget_check.json
"""
import json
import os
import re

import numpy as np
import pandas as pd

import common as K

TEX = os.path.join(K.HERE, "..", "paper", "tab_pooled.tex")


def _med(x):
    return float(np.median(x)) if len(x) else float("nan")


def table_s31():
    """The 2D-1 rows of Table S-3.1 as printed."""
    rows, on = [], False
    for line in open(TEX):
        if line.startswith("OpenFOAM 2D-1"):
            on = True
        elif on and not line.startswith(" &"):
            break
        if on:
            cells = [c.strip() for c in line.rstrip("\\\n").split("&")]
            nums = [c.replace("{,}", "") for c in cells[3:]]
            rows.append(dict(searches=int(nums[0]), unconv=int(nums[1]), accepted=int(nums[2]),
                             distinct=int(nums[3]), missed=int(nums[4]), declined=int(nums[5]),
                             pen_med=float(nums[7]), pen_max=float(nums[8]), runs=float(nums[9])))
    return rows


def main():
    out = {}
    printed = table_s31()
    stored = {(r["ob"], r["target"]): r for r in
              json.load(open(os.path.join(K.REPO, "paper", "data", "pooled_v1.2.1.json"))) if r["case"] == "2D-1"}
    base = json.load(open(os.path.join(K.REPO, "paper", "data", "baselines_v6.json")))["cells"]
    a1 = pd.read_csv(os.path.join(K.RESULTS, "arm1_rows.csv"))
    ok_all = True
    for (ob, tg), tab in zip(K.TESTS, printed):
        o = K.P.summarise(K.P.replay(K.CASE, ob, tg, K.CFG))
        mine = dict(searches=o["n"], unconv=o["start_failed"], accepted=o["accepted"], distinct=o["distinct"],
                    missed=o["missed"], declined=o["declined_reach"] + o["declined_unreach"],
                    pen_med=o["pen_med"], pen_max=o["pen_max"], runs=o["runs_med"])
        st = stored[(ob, tg)]
        g = a1[(a1.ob == ob) & (a1.target == tg)]
        met = g[g.eb_err <= 1]
        k = f"2D-1|{ob}|{tg}"
        t7 = dict(reached=int(len(met)), missed=int((g.eb_err > 1).sum()), declined=int(g.eb_q.isna().sum()),
                  pen=float(np.median(met.eb_cost / met.c_opt)), runs=float(np.median(g.runs_used)),
                  K=float(np.median(met.spent / met.eb_cost)), spent_median=float(np.median(g.spent)))
        sig3 = lambda a, b: f"{a:.3g}" == f"{b:.3g}"
        checks = dict(
            table_counts=all(mine[c] == tab[c] for c in ("searches", "unconv", "accepted", "distinct", "missed", "declined")),
            table_cost_ratio=sig3(mine["pen_med"], tab["pen_med"]) and sig3(mine["pen_max"], tab["pen_max"]),
            table_runs=mine["runs"] == tab["runs"],
            stored_json=all(mine[a] == st[b] for a, b in [("accepted", "accepted"), ("missed", "missed"),
                                                           ("distinct", "distinct"), ("runs", "runs_med")])
            and np.isclose(mine["pen_med"], st["pen_med"]) and np.isclose(mine["pen_max"], st["pen_max"]),
            exp72=t7["reached"] == base[k]["program"]["reached"] and np.isclose(t7["pen"], base[k]["program"]["pen"])
            and np.isclose(t7["K"], base[k]["program"]["K"]) and t7["runs"] == base[k]["program"]["runs"])
        ok_all &= all(checks.values())
        out[f"{ob} {tg:g}"] = dict(replay=mine, table_s31=tab, table7_replay=t7,
                                   table7_stored=base[k]["program"], checks=checks)
        print(f"{ob} {tg:<6g} S-3.1 {mine['searches']}/{mine['unconv']}/{mine['accepted']}/{mine['distinct']}/"
              f"{mine['missed']}/{mine['declined']} pen {mine['pen_med']:.3g}/{mine['pen_max']:.3g} runs {mine['runs']:.0f}"
              f" | printed {tab['searches']}/{tab['unconv']}/{tab['accepted']}/{tab['distinct']}/{tab['missed']}/"
              f"{tab['declined']} {tab['pen_med']:.3g}/{tab['pen_max']:.3g} {tab['runs']:.0f} | K {t7['K']:.1f} "
              f"(stored {base[k]['program']['K']:.1f}) | {checks}", flush=True)
    Ks = [v["table7_replay"]["K"] for v in out.values()]
    print(f"median of per-case K: {np.median(Ks):.0f} (Table 7: 224); all checks pass: {ok_all}")
    out["_summary"] = dict(median_K=float(np.median(Ks)), all_checks_pass=bool(ok_all))

    # errbudget-3p
    LAD3 = K.CASE.LAD[:3]
    full = lambda i: (i[0], i[1], i[2], K.NC_DEFAULT)
    st3 = [i[:3] for i in K.starts() if i[3] == K.NC_DEFAULT]
    g10 = [i for i in K.CASE.grid() if i[3] == K.NC_DEFAULT]
    res3 = {}
    for ob, tg in K.TESTS:
        c_opt = K.cost(K.optimum(ob, tg))
        c_opt10 = K.cost(K.optimum(ob, tg, grid=g10))
        acc = miss = dec = 0
        pen, pen10, Kk = [], [], []
        for a in st3:
            r, used = K.errbudget_search(ob, tg, a, ladders=LAD3, tensor_index=full)
            spent = sum(K.cost(i) for i in used if np.isfinite(K.cost(i)))
            if not r["certified"]:
                dec += 1
                continue
            q = full(tuple(r["proposal"]))
            acc += 1
            if K.grade(q, ob, tg) > 1:
                miss += 1
                continue
            pen.append(K.cost(q) / c_opt)
            pen10.append(K.cost(q) / c_opt10)
            Kk.append(spent / K.cost(q))
        res3[f"{ob} {tg:g}"] = dict(starts=len(st3), accepted=acc, missed=miss, declined=dec,
                                    pen_med=_med(pen), pen_max=float(np.max(pen)) if pen else float("nan"),
                                    pen10_med=_med(pen10), K_med=_med(Kk),
                                    c_opt_ratio_nc10=c_opt10 / c_opt)
        print(f"errbudget-3p {ob} {tg:<6g} starts {len(st3)} accepted {acc} missed {miss} declined {dec} "
              f"pen {_med(pen):.2f} (vs nc=10 optimum {_med(pen10):.2f}) K {_med(Kk):.0f}; "
              f"optimum at nc=10 / study optimum = {c_opt10 / c_opt:.3f}", flush=True)
    out["errbudget_3p"] = res3
    json.dump(out, open(os.path.join(K.RESULTS, "errbudget_check.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
