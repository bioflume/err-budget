"""Exploratory check of PROTOCOL.md section 7: does GRE change the outcome?

On 1 random Arm 1 start per test case (seed 1), the run set errbudget used is fitted with
  * SPRE, white-noise kernel (spre_white.py, the method of the comparison),
  * GRE as implemented in the authors' code (white-noise base kernel, stepwise rate),
and each posterior goes through the protocol's selector at z = 2.  (A Matern-1/2 SPRE fit was
dropped: one fit with the authors' code ran for over 10 minutes.)

    PYTHONPATH=ext/sparse-probabilistic-richardson-extrapolation/src python explore_kernels.py
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

import common as K
import spre_white as SW


def one(job):
    import jax.numpy as jnp
    from sparse_pre.extrapolation import extrapolation
    ob, tg, k = job
    a = K.starts()[k]
    r, used = K.errbudget_search(ob, tg, a)
    X, Y, _ = K.pre_dataset(used, ob)
    eps = K.eps_abs(ob, tg)
    out = dict(ob=ob, target=tg, start=list(a), n=len(Y),
               errbudget_err=K.grade(tuple(r["proposal"]), ob, tg) if r["certified"] else None)
    fits = {"SPRE white (fast)": lambda: SW.fit(X, Y)}
    for name, opts in [("GRE (authors' code)", {"name": "GRE", "k_name": "white"})]:
        def f(opts=opts):
            o = extrapolation(jnp.array(X), jnp.array(Y), dict(opts, plot=False))
            return dict(mu=float(np.ravel(o["mu"])[0]), sd=float(np.sqrt(max(float(np.ravel(o["var"])[0]), 0.0))))
        fits[name] = f
    for name, f in fits.items():
        t0 = time.time()
        try:
            res = f()
        except Exception as exc:                       # report, do not hide, a failed fit
            out[name] = dict(error=repr(exc))
            continue
        m, s = res["mu"], res["sd"]
        q = K.select(used, ob, tg, m, s, 2.0)
        out[name] = dict(m=m, s=s, est_err=abs(m - K.CASE.exact[ob]) / eps, s_over_eps=s / eps,
                         covered_2s=bool(abs(m - K.CASE.exact[ob]) <= 2 * s),
                         selected_err=K.grade(q, ob, tg) if q is not None else None,
                         selected_cost_ratio=K.cost(q) / K.cost(K.optimum(ob, tg)) if q is not None else None,
                         seconds=time.time() - t0)
    with open(os.path.join(K.RESULTS, "explore_kernels.jsonl"), "a") as fh:
        fh.write(json.dumps(out, default=float) + "\n")
    print(json.dumps({kk: (vv if not isinstance(vv, dict) else {a_: (round(b_, 4) if isinstance(b_, float) else b_)
                                                                 for a_, b_ in vv.items()}) for kk, vv in out.items()}),
          flush=True)
    return out


def main():
    rng = np.random.default_rng(1)
    n = len(K.starts())
    jobs = [(ob, tg, int(k)) for ob, tg in K.TESTS for k in rng.choice(n, 1, replace=False)]
    # Reduced after Arms 1-2 for the CPU limit: one test case per observable (REPORT.md, deviations)
    jobs = [j for j in jobs if (j[0], j[1]) in {("cd", 1e-2), ("dp", 1e-2), ("cl", 3e-2)}]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        rows = pool.map(one, jobs)
    json.dump(rows, open(os.path.join(K.RESULTS, "explore_kernels.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
