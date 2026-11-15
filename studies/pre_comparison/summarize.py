"""Outcome measures of PROTOCOL.md section 5 and the headline criteria of section 6.

    python summarize.py      # reads results/arm1_rows.csv, results/arm2_rows.csv,
                             # results/errbudget_check.json; writes results/summary.json
                             # and results/tables.md
"""
import json
import os

import numpy as np
import pandas as pd

import common as K

LAB = {"cd": "c_D", "dp": "Δp", "cl": "c_L"}


def med(x):
    x = np.asarray(x, dtype=float)
    return float(np.median(x)) if len(x) else float("nan")


def pre_stats(g, z, spent_col="spent"):
    key = f"z{z:g}"
    err, c = g[f"pre_err_{key}"], g[f"pre_cost_{key}"]
    ret = err.notna()
    met = ret & (err <= 1)
    eps = np.array([K.eps_abs(ob, tg) for ob, tg in zip(g.ob, g.target)])
    return dict(n=int(len(g)), returned=int(ret.sum()), missed=int((err > 1).sum()),
                declined=int((~ret).sum()), worst_err=float(err.max()) if ret.any() else float("nan"),
                pen_med=med(c[met] / g.c_opt[met]), pen_max=float((c[met] / g.c_opt[met]).max()) if met.any() else float("nan"),
                pen_med_all=med(c[ret] / g.c_opt[ret]), K_med=med(g[spent_col][met] / c[met]))


def estimate_stats(g):
    eps = np.array([K.eps_abs(ob, tg) for ob, tg in zip(g.ob, g.target)])
    ok = g.pre_m.notna() & g.pre_s.notna()
    e = g.pre_est_err[ok].to_numpy()
    s_eps = (g.pre_s[ok] / eps[ok.to_numpy()]).to_numpy()
    return dict(fits=int(ok.sum()), est_err_med=med(e), est_err_p90=float(np.percentile(e, 90)) if len(e) else np.nan,
                est_err_max=float(e.max()) if len(e) else np.nan,
                coverage={f"z{z:g}": float(np.mean(e <= z * s_eps)) for z in K.ZS},
                frac_2s_le_eps=float(np.mean(2 * s_eps <= 1)), frac_m_within_eps=float(np.mean(e <= 1)),
                s_over_eps_med=med(s_eps))


def eb_stats(g):
    met = g.eb_err <= 1
    return dict(n=int(len(g)), missed=int((g.eb_err > 1).sum()), declined=int(g.eb_q.isna().sum()),
                pen_med=med(g.eb_cost[met] / g.c_opt[met]), pen_max=float((g.eb_cost[met] / g.c_opt[met]).max()),
                K_med=med(g.spent[met] / g.eb_cost[met]), runs_med=med(g.runs_used))


def main():
    a1 = pd.read_csv(os.path.join(K.RESULTS, "arm1_rows.csv"), keep_default_na=True)
    a2 = pd.read_csv(os.path.join(K.RESULTS, "arm2_rows.csv"), keep_default_na=True)
    chk = json.load(open(os.path.join(K.RESULTS, "errbudget_check.json")))
    S = dict(arm1={}, arm2={})
    for ob, tg in K.TESTS:
        key = f"{ob} {tg:g}"
        g = a1[(a1.ob == ob) & (a1.target == tg)]
        S["arm1"][key] = dict(errbudget=eb_stats(g), pre={f"z{z:g}": pre_stats(g, z) for z in K.ZS},
                              estimate=estimate_stats(g),
                              pre_n_med=med(g.pre_n), A_mode=g.pre_A.mode().iloc[0])
        for fac in sorted(a2.budget_factor.unique()):
            h = a2[(a2.ob == ob) & (a2.target == tg) & (a2.budget_factor == fac)]
            S["arm2"].setdefault(f"B x{fac:g}", {})[key] = dict(
                budget=float(h.budget.iloc[0]), spent_med=med(h.spent), runs_med=med(h.runs),
                converged_runs_med=med(h.converged_runs),
                pre={f"z{z:g}": pre_stats(h, z) for z in K.ZS}, estimate=estimate_stats(h),
                errbudget=eb_stats(g), errbudget_3p=chk["errbudget_3p"][key])
    S["headline"] = headline(S)
    json.dump(S, open(os.path.join(K.RESULTS, "summary.json"), "w"), indent=1, default=float)
    write_tables(S)


def headline(S):
    """Section 6 of the protocol, at z = 2."""
    out = {}
    for arm, rows in [("arm1", S["arm1"]), ("arm2", S["arm2"]["B x1"])]:
        miss_where_eb_not = [k for k, v in rows.items()
                             if v["pre"]["z2"]["missed"] >= max(1, 0.01 * v["pre"]["z2"]["n"])
                             and v["errbudget"]["missed"] == 0]
        low_cov = [k for k, v in rows.items() if v["estimate"]["coverage"]["z2"] < 0.90]
        no_miss = all(v["pre"]["z2"]["missed"] == 0 for v in rows.values())
        decl_ok = all(v["pre"]["z2"]["declined"] <= 0.05 * v["pre"]["z2"]["n"] for v in rows.values())
        cheaper = [k for k, v in rows.items() if v["pre"]["z2"]["pen_med"] * 1.5 <= v["errbudget"]["pen_med"]]
        many_decl = [k for k, v in rows.items() if v["pre"]["z2"]["declined"] > 0.10 * v["pre"]["z2"]["n"]]
        c1 = bool(miss_where_eb_not) or len(low_cov) >= 2
        c2 = no_miss and decl_ok and len(cheaper) >= 5
        if arm == "arm2":
            c2 = c2 and all(v["pre"]["z2"]["K_med"] <= v["errbudget"]["K_med"] for v in rows.values())
        c3 = no_miss and len(many_decl) >= 3
        out[arm] = dict(criterion1_less_reliable=c1, pre_misses_where_errbudget_does_not=miss_where_eb_not,
                        coverage_below_90=low_cov, criterion2_cheaper_same_reliability=c2,
                        cases_cost_ratio_1p5x_lower=cheaper, criterion3_less_usable=c3,
                        cases_declining_over_10pct=many_decl)
    return out


def fmt(x, d=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    if abs(x) >= 100:
        return f"{x:,.0f}"
    return f"{x:.{d}g}" if abs(x) < 10 else f"{x:.3g}"


def write_tables(S):
    L = []
    L.append("## Arm 1: PRE sees exactly errbudget's runs (608 searches per case)\n")
    L.append("| Test case | errbudget missed / declined | PRE missed / declined | errbudget cost ratio | PRE cost ratio | "
             "errbudget K | PRE K | PRE coverage of Q* (±2s) | median |m−Q*|/ε |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for k, v in S["arm1"].items():
        ob, tg = k.split()
        e, p, est = v["errbudget"], v["pre"]["z2"], v["estimate"]
        L.append(f"| {LAB[ob]} {tg} | {e['missed']} / {e['declined']} | {p['missed']} / {p['declined']} | "
                 f"{fmt(e['pen_med'])} | {fmt(p['pen_med'])} | {fmt(e['K_med'])} | {fmt(p['K_med'])} | "
                 f"{100 * est['coverage']['z2']:.0f}% | {fmt(est['est_err_med'])} |")
    L.append("\nPRE at z = 1 and z = 3 (missed / declined, median cost ratio):\n")
    L.append("| Test case | z = 1 | z = 3 | coverage ±1s / ±3s | 2s ≤ ε | |m−Q*| ≤ ε |")
    L.append("|---|---|---|---|---|---|")
    for k, v in S["arm1"].items():
        ob, tg = k.split()
        p1, p3, est = v["pre"]["z1"], v["pre"]["z3"], v["estimate"]
        L.append(f"| {LAB[ob]} {tg} | {p1['missed']} / {p1['declined']}, {fmt(p1['pen_med'])} | "
                 f"{p3['missed']} / {p3['declined']}, {fmt(p3['pen_med'])} | {100 * est['coverage']['z1']:.0f}% / "
                 f"{100 * est['coverage']['z3']:.0f}% | {100 * est['frac_2s_le_eps']:.0f}% | {100 * est['frac_m_within_eps']:.0f}% |")
    for fac, rows in S["arm2"].items():
        L.append(f"\n## Arm 2: PRE's own design, budget {fac} errbudget's median search cost (50 seeds per case)\n")
        L.append("| Test case | budget (cells×iter.) | errbudget missed / declined | PRE missed / declined | errbudget cost ratio | "
                 "PRE cost ratio | errbudget-3p cost ratio | errbudget K | PRE K | PRE coverage (±2s) | median |m−Q*|/ε |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for k, v in rows.items():
            ob, tg = k.split()
            e, p, est, e3 = v["errbudget"], v["pre"]["z2"], v["estimate"], v["errbudget_3p"]
            L.append(f"| {LAB[ob]} {tg} | {v['budget']:.3g} | {e['missed']}/608 / {e['declined']} | {p['missed']}/50 / {p['declined']} | "
                     f"{fmt(e['pen_med'])} | {fmt(p['pen_med'])} | {fmt(e3['pen_med'])} | {fmt(e['K_med'])} | {fmt(p['K_med'])} | "
                     f"{100 * est['coverage']['z2']:.0f}% | {fmt(est['est_err_med'])} |")
    L.append("\n## Headline criteria (protocol section 6, z = 2)\n")
    L.append("```\n" + json.dumps(S["headline"], indent=1) + "\n```")
    open(os.path.join(K.RESULTS, "tables.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
