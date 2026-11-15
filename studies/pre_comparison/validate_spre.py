"""Validation of PRE before it is used on the 2D-1 study.

    python validate_spre.py published   # step V1: the authors' code reproduces the published
                                        #   cubature example (Figs. 1-2 of arXiv 2604.02072)
    python validate_spre.py agree       # step V2: spre_white.py reproduces the authors' code
                                        #   on 2D-1 run sets (no grading against exact values)

V1 compares a fresh run of the authors' run_cubature_experiments.py (validation/rerun/) with the
results stored in their repository, from which the paper's Figs. 1 and 2 are drawn, and checks
the convergence rates the paper states for SPRE with the white-noise kernel (Section 3.2, after
Proposition 1): O(h^2) for s = 0 and O(h^4) for s = 1, until the floating-point floor.

V2 fits the same data with the authors' extrapolation() (SPRE, white-noise kernel, stepwise
selection of A, trust-region optimisation of the amplitude) and with spre_white.fit, on 21 run
sets of errbudget's searches (3 random starts per test case, seed 0), and on the n = 4 SPRE
data set of Fig. 1.  Writes results/validation.json.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(HERE, "ext", "sparse-probabilistic-richardson-extrapolation")
PUB = os.path.join(EXT, "data", "cubature", "results")
NEW = os.path.join(HERE, "validation", "rerun", "data", "cubature", "results")
OUT = os.path.join(HERE, "results", "validation.json")


def _load(path):
    return json.load(open(path)) if os.path.exists(path) else {}


def _save(key, val):
    doc = _load(OUT)
    doc[key] = val
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(doc, open(OUT, "w"), indent=1, default=float)


def published():
    import pandas as pd
    pub = pd.read_csv(os.path.join(PUB, "convergence.dat"), sep="\t")
    new = pd.read_csv(os.path.join(NEW, "convergence.dat"), sep="\t")
    assert len(pub) == len(new), (len(pub), len(new))
    out = {}
    # 1. the rerun against the stored (published) numbers
    # Rows whose error is at the floating-point floor (< 1e-12) differ by round-off only; they
    # are compared in absolute terms.  GRE (Matern-3/2, 10 trust-region iterations) is not used
    # in the comparison and is reported for completeness.
    for col in ["abs_err_baseline", "abs_err_MRE", "abs_err_SPRE", "rel_err_SPRE",
                "abs_err_GRE", "rel_err_GRE"]:
        a, b = pub[col].to_numpy(), new[col].to_numpy()
        base = pub[col.replace("rel_err", "abs_err")].to_numpy()
        above = base > 1e-12
        rel = np.abs(a - b) / np.maximum(np.abs(a), 1e-300)
        out[col] = dict(rows=len(rel), rows_above_floor=int(above.sum()),
                        max_rel_diff_above_floor=float(rel[above].max()),
                        max_abs_diff_at_floor=float(np.abs(a - b)[~above].max()) if (~above).any() else 0.0)
        print(f"V1 {col:<18} rerun vs published, {above.sum()} rows above the floor: max rel. diff "
              f"{rel[above].max():.1e}; at the floor max abs. diff {out[col]['max_abs_diff_at_floor']:.1e}")
    # 2. the rates stated in the paper, measured on the rerun: least-squares slope of
    #    log(abs error) against log(h) over the rows above the floating-point floor
    rates = {}
    for (d, s), g in new.groupby(["d", "s"]):
        g = g[(g["abs_err_SPRE"] > 1e-13) & (g["h"] <= 0.5)]
        if len(g) >= 3:
            slope = np.polyfit(np.log(g["h"]), np.log(g["abs_err_SPRE"]), 1)[0]
        else:
            slope = np.nan
        expect = 2 if s == 0 else 4
        rates[f"d={d},s={s}"] = dict(slope=float(slope), stated=expect, rows=int(len(g)))
        print(f"V1 rate SPRE d={d} s={s}: measured {slope:.2f} over {len(g)} rows; paper: O(h^{expect})")
    # 3. calibration claim: actual error usually within 3 standard deviations
    within3 = float(np.mean(np.abs(new["rel_err_SPRE"]) <= 3))
    print(f"V1 SPRE |error| / sd <= 3 in {100 * within3:.0f}% of rows (paper: 'usually')")
    # 4. Fig. 1 (illustration): estimates of the three methods
    ip, inew = _load(os.path.join(PUB, "illustration.json")), _load(os.path.join(NEW, "illustration.json"))
    ill = {}
    for m in ["MRE", "GRE", "SPRE"]:
        ill[m] = dict(published=ip[m]["mu"], rerun=inew[m]["mu"],
                      rel_diff=abs(ip[m]["mu"] - inew[m]["mu"]) / abs(ip[m]["mu"]))
        print(f"V1 Fig. 1 {m}: published {ip[m]['mu']:.10f}, rerun {inew[m]['mu']:.10f} "
              f"(truth {ip['f0']:.10f})")
    _save("V1_published", dict(columns=out, rates=rates, within3=within3, illustration=ill))


def agree():
    """Run 1 (seed 0, 3 starts per case, in the order of common.TESTS) was stopped after the
    c_D 3e-3 case for CPU time (the authors' code took up to 10 min per 2D-1 run set); its five
    rows are in validation/agree_run1.log.  Run 2 (`agree run2`) takes 1 start per remaining
    case (seed 2).  Rows are saved after each data set."""
    import jax.numpy as jnp
    from sparse_pre.extrapolation import extrapolation
    import common as K
    import spre_white as SW
    run2 = len(sys.argv) > 2 and sys.argv[2] == "run2"
    rng = np.random.default_rng(2 if run2 else 0)
    st = K.starts()
    rows = []
    key = "V2_agreement_run2" if run2 else "V2_agreement"

    def compare(label, X, Y):
        t0 = time.time()
        o = extrapolation(jnp.array(X), jnp.array(Y), {"name": "SPRE", "k_name": "white", "plot": False})
        ta = time.time() - t0
        mu_a = float(np.ravel(o["mu"])[0])
        sd_a = float(np.sqrt(np.ravel(o["var"])[0]))
        A_a = sorted(map(tuple, np.asarray(o["A"]).tolist()))
        f = SW.fit(X, Y)
        A_f = sorted(map(tuple, f["A"].tolist()))
        r = dict(label=label, n=len(Y), mu_authors=mu_a, mu_fast=f["mu"], sd_authors=sd_a,
                 sd_fast=f["sd"], A_authors=A_a, A_fast=A_f, same_A=A_a == A_f,
                 rel_diff_mu_over_sd=abs(mu_a - f["mu"]) / sd_a, rel_diff_sd=abs(sd_a - f["sd"]) / sd_a,
                 seconds_authors=ta)
        rows.append(r)
        _save(key, rows)
        print(f"V2 {label:<24} n={len(Y):<4} same A {r['same_A']}  |dmu|/sd {r['rel_diff_mu_over_sd']:.1e}  "
              f"|dsd|/sd {r['rel_diff_sd']:.1e}  ({ta:.0f} s)", flush=True)

    if run2:
        for ob, tg in K.TESTS[2:]:
            k = int(rng.choice(len(st)))
            r, used = K.errbudget_search(ob, tg, st[k])
            X, Y, _ = K.pre_dataset(used, ob)
            compare(f"2D-1 {ob} {tg:g} start {k}", X, Y)
        return
    ill = _load(os.path.join(PUB, "illustration.json"))["SPRE"]
    compare("Fig. 1 SPRE data", np.array(ill["X"]), np.array(ill["Y"]))
    for ob, tg in K.TESTS:
        for k in rng.choice(len(st), 3, replace=False):
            r, used = K.errbudget_search(ob, tg, st[k])
            X, Y, _ = K.pre_dataset(used, ob)
            compare(f"2D-1 {ob} {tg:g} start {k}", X, Y)
    print(f"V2 summary: same A in {sum(r['same_A'] for r in rows)}/{len(rows)}; "
          f"max |dmu|/sd {max(r['rel_diff_mu_over_sd'] for r in rows):.1e}; "
          f"max |dsd|/sd {max(r['rel_diff_sd'] for r in rows):.1e}")


if __name__ == "__main__":
    {"published": published, "agree": agree}[sys.argv[1]]()
