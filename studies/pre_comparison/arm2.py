"""Arm 2 of PROTOCOL.md: PRE chooses its own runs under errbudget's median search cost.

The budget B of a test case is the median, over errbudget's 608 searches, of the summed cost of
the runs it used (from results/arm1_rows.csv).  PRE spends it with the sequential cost-constrained
design of the sparse PRE paper (Section 3.3, eq. 11; the authors' run_sequential_design.py), on the
210 settings at the default nCellsInCoarsestLevel:

  1. an 8-point Latin hypercube over the level indices; its 4 cheapest settings are run;
  2. three iterations, each: fit SPRE to the converged runs; draw 1,000 random batches within the
     iteration's budget; run the batch with the smallest posterior variance of Q*;
     unspent budget carries over;
  3. final SPRE fit and the protocol's selector over the converged runs of the design.

Costs are taken as known in advance (stored costs).  50 seeds per test case, at B and 3B.

    python arm2.py            # writes results/arm2_rows.csv
"""
import csv
import os
import time
from multiprocessing import Pool

import numpy as np

import common as K
import spre_white as SW

N_SEEDS = 50
N_ITER = 3
N_BATCHES = 1000
BUDGET_FACTORS = (1.0, 3.0)
CAND = [i for i in K.CASE.grid() if i[3] == K.NC_DEFAULT and K.usable(i)]
LEVELS = [len(l) for l in K.CASE.LAD[:3]]


def latin_hypercube(rng, n=8):
    """n distinct settings: each axis split into n strata, one uniform draw per stratum,
    mapped to the axis's level indices."""
    for _ in range(1000):
        cols = []
        for m in LEVELS:
            u = (rng.permutation(n) + rng.uniform(size=n)) / n
            cols.append(np.minimum((u * m).astype(int), m - 1))
        pts = {(int(a), int(b), int(c), K.NC_DEFAULT) for a, b, c in zip(*cols)}
        if len(pts) == n:
            return sorted(pts)
    raise RuntimeError("no distinct Latin hypercube")


def norm_inputs(X_ref):
    """The normalisation SPRE applies to the current data (range of each input)."""
    if len(X_ref) >= 2:
        return SW.EP + (X_ref.max(axis=0) - X_ref.min(axis=0))
    return np.abs(np.asarray([K.x_of(i) for i in CAND])).max(axis=0)


def design(ob, tg, seed, B):
    rng = np.random.default_rng(seed)
    runs = sorted(latin_hypercube(rng), key=K.cost)[:4]
    spent = sum(K.cost(i) for i in runs)
    per = max(B - spent, 0.0) / N_ITER
    carry = 0.0
    for it in range(N_ITER):
        bud = per + carry
        X, Y, _ = K.pre_dataset(runs, ob)
        f = SW.fit(X, Y) if len(Y) >= 2 else None
        A = f["A"] if (f is not None and np.isfinite(f["mu"])) else np.zeros((1, 3), int)
        nX = norm_inputs(X)
        unused = [i for i in CAND if i not in runs]
        best, best_val, best_cost = None, np.inf, 0.0
        fallback = None
        for _ in range(N_BATCHES):
            batch, s = [], 0.0
            for j in rng.permutation(len(unused)):
                c = K.cost(unused[j])
                if s + c <= bud:
                    batch.append(unused[j])
                    s += c
                elif batch:
                    break
            if not batch:
                continue
            Xall = np.vstack([X.reshape(-1, 3), np.array([K.x_of(i) for i in batch])]) / nX
            val = SW.design_variance(Xall, A)
            if val < best_val:
                best, best_val, best_cost = batch, val, s
            if fallback is None or len(batch) > len(fallback[0]):
                fallback = (batch, s)
        if best is None and fallback is not None:
            best, best_cost = fallback
        if best is None:
            carry = bud
            continue
        runs += best
        spent += best_cost
        carry = bud - best_cost
    X, Y, _ = K.pre_dataset(runs, ob)
    f = SW.fit(X, Y) if len(Y) >= 2 else dict(mu=np.nan, sd=np.nan, A=np.zeros((1, 3), int))
    row = dict(ob=ob, target=tg, seed=seed, budget=B, spent=spent, runs=len(runs),
               converged_runs=len(Y), pre_m=f["mu"], pre_s=f["sd"],
               pre_A=";".join("".join(map(str, a)) for a in np.asarray(f["A"]).tolist()),
               pre_est_err=abs(f["mu"] - K.CASE.exact[ob]) / K.eps_abs(ob, tg),
               c_opt=K.cost(K.optimum(ob, tg)), finest_run=max(K.x_of(i)[0] ** -1 for i in runs))
    for z in K.ZS:
        q = K.select(runs, ob, tg, f["mu"], f["sd"], z)
        key = f"z{z:g}"
        if q is None:
            row.update({f"pre_q_{key}": "", f"pre_err_{key}": np.nan, f"pre_cost_{key}": np.nan})
        else:
            row.update({f"pre_q_{key}": "-".join(map(str, q)), f"pre_err_{key}": K.grade(q, ob, tg),
                        f"pre_cost_{key}": K.cost(q)})
    return row


def run_case(args):
    (ob, tg), B0 = args
    t0 = time.time()
    rows = []
    for fac in BUDGET_FACTORS:
        for seed in range(N_SEEDS):
            r = design(ob, tg, seed, fac * B0)
            r["budget_factor"] = fac
            rows.append(r)
    print(f"arm 2 {ob} {tg:g}: B = {B0:.3g}, {len(rows)} designs in {time.time() - t0:.0f} s", flush=True)
    return rows


def budgets():
    """Median summed cost of errbudget's searches per test case, from Arm 1."""
    import pandas as pd
    a1 = pd.read_csv(os.path.join(K.RESULTS, "arm1_rows.csv"))
    return {(ob, tg): float(np.median(g["spent"])) for (ob, tg), g in a1.groupby(["ob", "target"])}


def main():
    Bs = budgets()
    jobs = [((ob, tg), Bs[(ob, tg)]) for ob, tg in K.TESTS]
    with Pool(len(jobs)) as pool:
        parts = pool.map(run_case, jobs)
    rows = [r for p in parts for r in p]
    path = os.path.join(K.RESULTS, "arm2_rows.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path}: {len(rows)} rows")


if __name__ == "__main__":
    main()
