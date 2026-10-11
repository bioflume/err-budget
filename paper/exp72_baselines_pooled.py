"""The program and three procedures in common use, graded as the pooled table
(referee points 15, 16 and 20).

For every admitted cell of the selected cases and every start whose own run
converged, each procedure returns a setting.  A returned setting *reaches* the
target if its run converged and its true error is within the target; otherwise it
*misses* (the program may also decline, which is neither).  Runs whose nested
solver did not converge are failed runs for every procedure.  Reported per case and
pooled:
  reached, missed, declined;
  penalty: cost of the returned setting / cheapest converged setting that reaches
           the target (median over the cells of the per-cell median);
  runs:    distinct solver evaluations used to decide (median of cell medians);
  K:       calibration cost in production runs, i.e. the summed cost of the runs used
           to decide divided by the cost of one run at the returned setting
           (median of cell medians) -- the same unit for every procedure.
    python exp72_baselines_pooled.py [CASE ...]
    KNOBS_BASE_OUT=baselines_lmp_gmx.json python exp72_baselines_pooled.py LMP-20 LMP-40 LMP-64 GMX-d
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C
import pooled as P
from baselines import seq, coord, ref_check
from knobs import run_protocol


def counting(look):
    seen = {}

    def f(i):
        i = tuple(i)
        if i not in seen:
            seen[i] = look(i)
        return seen[i]
    f.seen = seen
    return f


def cell(case, ob, tg, cfg):
    eps = tg * abs(case.exact[ob])
    T = case.tensor
    G = case.grid()
    good = [i for i in G if case.usable(i) and case.converged(i)]
    err = {i: abs(T[i][0][ob] - case.exact[ob]) for i in good}
    feas = [i for i in good if err[i] <= eps]
    if not feas:
        return None
    c_opt = min(T[i][1] for i in feas)

    def look(i):
        i = tuple(i)
        if not case.usable(i) or not case.converged(i):
            return np.nan, (T[i][1] if np.isfinite(T[i][1]) else np.inf)
        return T[i][0][ob], T[i][1]

    procs = {"program": None, "seq-1": 1.0, "seq-3": 3.0, "seq-10": 10.0, "coord-10": 10.0,
             "ref": None}
    out = {p: dict(reached=0, missed=0, declined=0, pen=[], runs=[], K=[]) for p in procs}
    kw = dict(cfg["kw"])
    for a in good:
        for p, div in procs.items():
            lk = counting(look)
            if p == "program":
                r = run_protocol(lk, case.LAD, a, eps, hs=case.HS, kinds=case.KINDS,
                                 models=case.models, **kw)
                if not r["certified"]:
                    out[p]["declined"] += 1
                    out[p]["runs"].append(len(lk.seen))
                    continue
                q = tuple(r["proposal"])
            elif p.startswith("seq"):
                q, _ = seq(lk, case.LAD, a, eps / div)
            elif p.startswith("coord"):
                q, _ = coord(lk, case.LAD, a, eps / div)
            else:
                q, _ = ref_check(lk, case.LAD, good, eps)
                for i in good:
                    lk(i)
            q = tuple(q)
            o = out[p]
            o["runs"].append(len(lk.seen))
            spent = sum(v[1] for v in lk.seen.values() if np.isfinite(v[1]))
            if q in err and err[q] <= eps:
                o["reached"] += 1
                o["pen"].append(T[q][1] / c_opt)
                o["K"].append(spent / T[q][1])
            else:
                o["missed"] += 1
    return out


def main():
    keys = sys.argv[1:] or ["L2", "L3", "L4", "L5", "2D-1", "2D-2", "QE"]
    cfg = P.CONFIGS[os.environ.get("KNOBS_CFG", "v1.2.1")]
    loaders = [C._manufactured, lambda: [C._dfg1()], lambda: [C._dfg2()], lambda: [C._qe()],
               C._lammps, C._gromacs, lambda: [C._abinit()]]
    cases = [c for L in loaders for c in L() if c.key in keys]
    res = {}
    for case in cases:
        for ob, tgs in case.targets.items():
            for tg in tgs:
                if not case.admitted(ob, tg):
                    continue
                o = cell(case, ob, tg, cfg)
                if o is None:
                    continue
                res[f"{case.key}|{ob}|{tg}"] = {p: dict(reached=v["reached"], missed=v["missed"],
                                                        declined=v["declined"],
                                                        pen=float(np.median(v["pen"])) if v["pen"] else np.nan,
                                                        runs=float(np.median(v["runs"])),
                                                        K=float(np.median(v["K"])) if v["K"] else np.nan)
                                                for p, v in o.items()}
                print(case.key, ob, tg, {p: (v["reached"], v["missed"], v["declined"],
                                             round(float(np.median(v["pen"])), 2) if v["pen"] else None,
                                             int(np.median(v["runs"])))
                                         for p, v in o.items()}, flush=True)
    out = "baselines_pooled.json" if not sys.argv[1:] else os.environ.get("KNOBS_BASE_OUT", "baselines_pooled.json")
    json.dump(res, open(os.path.join(C.D, out), "w"), indent=1)
    # pooled per case group
    groups = {"L2-L4": ["L2", "L3", "L4"], "L5": ["L5"], "2D-1": ["2D-1"], "2D-2": ["2D-2"], "QE": ["QE"]}
    for g, ks in groups.items():
        rows = {k: v for k, v in res.items() if k.split("|")[0] in ks}
        if not rows:
            continue
        print(f"\n{g}: {len(rows)} cells")
        for p in ["program", "seq-1", "seq-3", "seq-10", "coord-10", "ref"]:
            R = [v[p] for v in rows.values()]
            pen = [r["pen"] for r in R if np.isfinite(r["pen"])]
            Ks = [r["K"] for r in R if np.isfinite(r["K"])]
            print(f"  {p:<9} reached {sum(r['reached'] for r in R):>6} missed {sum(r['missed'] for r in R):>5} "
                  f"declined {sum(r['declined'] for r in R):>5} penalty {np.median(pen) if pen else np.nan:6.2f} "
                  f"runs {np.median([r['runs'] for r in R]):5.0f} K {np.median(Ks) if Ks else np.nan:7.1f}")


if __name__ == "__main__":
    main()
