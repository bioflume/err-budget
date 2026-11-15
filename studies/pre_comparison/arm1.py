"""Arm 1 of PROTOCOL.md: PRE sees exactly the runs errbudget's search used.

For each graded test case and each of the 608 converged starts, errbudget (v1.2.1) is replayed
as in exp72_baselines_pooled.cell and every run it asks for is recorded.  SPRE (white-noise
kernel, spre_white.py) is fitted to the converged runs among them; the selector of the protocol
returns the cheapest seen converged run with |Q - m| + z s <= eps, for z = 1, 2, 3.

    python arm1.py            # writes results/arm1_rows.csv (one row per test case and start)
"""
import csv
import os
import time
from multiprocessing import Pool

import numpy as np

import common as K
import spre_white as SW


def run_case(test):
    ob, tg = test
    t0 = time.time()
    opt = K.optimum(ob, tg)
    c_opt = K.cost(opt)
    rows = []
    for a in K.starts():
        r, used = K.errbudget_search(ob, tg, a)
        spent = sum(K.cost(i) for i in used if np.isfinite(K.cost(i)))
        row = dict(ob=ob, target=tg, start="-".join(map(str, a)), c_opt=c_opt,
                   runs_used=len(used), spent=spent)
        if r["certified"]:
            q = tuple(r["proposal"])
            row.update(eb_q="-".join(map(str, q)), eb_err=K.grade(q, ob, tg), eb_cost=K.cost(q))
        else:
            row.update(eb_q="", eb_err=np.nan, eb_cost=np.nan)
        X, Y, _ = K.pre_dataset(used, ob)
        row["pre_n"] = len(Y)
        f = SW.fit(X, Y) if len(Y) >= 2 else dict(mu=np.nan, sd=np.nan, A=np.zeros((1, 3), int))
        row.update(pre_m=f["mu"], pre_s=f["sd"],
                   pre_A=";".join("".join(map(str, a_)) for a_ in np.asarray(f["A"]).tolist()),
                   pre_est_err=abs(f["mu"] - K.CASE.exact[ob]) / K.eps_abs(ob, tg))
        for z in K.ZS:
            q = K.select(used, ob, tg, f["mu"], f["sd"], z)
            key = f"z{z:g}"
            if q is None:
                row.update({f"pre_q_{key}": "", f"pre_err_{key}": np.nan, f"pre_cost_{key}": np.nan})
            else:
                row.update({f"pre_q_{key}": "-".join(map(str, q)), f"pre_err_{key}": K.grade(q, ob, tg),
                            f"pre_cost_{key}": K.cost(q)})
        row["pre_degree"] = int(np.asarray(f["A"]).sum(axis=1).max())
        # sensitivity of the protocol: stepwise selection capped at order 4
        f4 = SW.fit(X, Y, max_order=4) if len(Y) >= 2 else dict(mu=np.nan, sd=np.nan)
        q = K.select(used, ob, tg, f4["mu"], f4["sd"], 2.0)
        row.update(pre4_m=f4["mu"], pre4_s=f4["sd"],
                   pre4_est_err=abs(f4["mu"] - K.CASE.exact[ob]) / K.eps_abs(ob, tg),
                   pre4_q_z2="-".join(map(str, q)) if q is not None else "",
                   pre4_err_z2=K.grade(q, ob, tg) if q is not None else np.nan,
                   pre4_cost_z2=K.cost(q) if q is not None else np.nan)
        rows.append(row)
    print(f"arm 1 {ob} {tg:g}: {len(rows)} starts in {time.time() - t0:.0f} s", flush=True)
    return rows


def main():
    os.makedirs(K.RESULTS, exist_ok=True)
    with Pool(len(K.TESTS)) as pool:
        parts = pool.map(run_case, K.TESTS)
    rows = [r for p in parts for r in p]
    path = os.path.join(K.RESULTS, "arm1_rows.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path}: {len(rows)} rows")


if __name__ == "__main__":
    main()
