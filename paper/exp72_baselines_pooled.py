"""The program and the three procedures of current practice (baselines.py) from every
converged starting setting of every graded test case (Section 6.3, Fig. 4, Tables S-4.2
and S-4.3).

    ERRBUDGET_CFG=v1.2.1 python exp72_baselines_pooled.py                  # every case
    ERRBUDGET_CFG=v1.2.1 ERRBUDGET_BASE_OUT=baselines_lmp_gmx.json \
        python exp72_baselines_pooled.py LMP-20 LMP-40 LMP-64 GMX-d       # Table S-4.3

A returned setting meets the target if its run converged and its true error is within the
target; otherwise it misses (the program may also decline).  Per test case and procedure:
met, missed, declined; the cost ratio against the optimal setting; the runs used to decide;
and the search cost K, the summed cost of those runs over the cost of one run at the
returned setting.  S and SJ take the parameters in the order of the case's list
(cases.py), discretization parameters first; GROMACS: rcoulomb, fourierspacing, pme-order,
ewald-rtol.  Writes data/baselines_v6.json (or ERRBUDGET_BASE_OUT).
"""
import json
import os
import sys

import numpy as np

import cases as C
import pooled as P
from baselines import seq, coord, ref_check
from errbudget import run_protocol

PROCS = {"program": None, "seq-1": 1.0, "seq-3": 3.0, "seq-10": 10.0,
         "coord-1": 1.0, "coord-3": 3.0, "coord-10": 10.0, "ref": None}
GROUPS = [("L2-L4", ["L2", "L3", "L4"]), ("L5", ["L5"]), ("2D-1", ["2D-1"]),
          ("2D-2 (rev. 3)", ["2D-2"]), ("2D-2 (refined)", ["2D-2h"]),
          ("LAMMPS", ["LMP-20", "LMP-40", "LMP-64"]), ("GROMACS double", ["GMX-d"]),
          ("GROMACS mixed", ["GMX-s"]), ("ABINIT", ["ABINIT"]), ("QE", ["QE"])]


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

    out = {p: dict(reached=0, missed=0, declined=0, pen=[], runs=[], K=[]) for p in PROCS}
    kw = dict(cfg["kw"])
    if case.floor is not None and not cfg.get("no_floor"):
        kw.setdefault("noise_floor", case.floor)       # as in Table 3 (pooled.replay)
    kinds = [cfg.get("kinds_map", {}).get(k, k) for k in case.KINDS]
    for a in good:
        for p, div in PROCS.items():
            lk = counting(look)
            if p == "program":
                r = run_protocol(lk, case.LAD, a, eps, hs=case.HS, kinds=kinds,
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


def load_cases(keys):
    loaders = [C._manufactured, lambda: [C._dfg1()], lambda: [C._dfg2()], lambda: [C._qe()],
               C._lammps, C._gromacs, lambda: [C._abinit()], lambda: [C._dfg2_refined()]]
    out = []
    for L in loaders:
        try:
            out += [c for c in L() if keys is None or c.key in keys]
        except FileNotFoundError as exc:
            print(f"skipping: {exc}", flush=True)
    return out


def main():
    keys = sys.argv[1:] or None
    cfgname = os.environ.get("ERRBUDGET_CFG", "v1.2.1")
    cfg = P.CONFIGS[cfgname]
    res = {}
    for case in load_cases(keys):
        print(f"{case.key}: sweep order of the baselines = parameter order "
              f"{[k[0] for k in case.KINDS]} (d discretisation, n nested solver, c cost-only)",
              flush=True)
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
    # one file for every case: a run on selected cases replaces only their cells
    path = os.path.join(C.D, os.environ.get("ERRBUDGET_BASE_OUT", "baselines_v6.json"))
    doc = dict(config=cfgname, cells={})
    if os.path.exists(path):
        old = json.load(open(path))
        if old.get("config") == cfgname:
            doc = old
    doc["cells"].update(res)
    json.dump(doc, open(path, "w"), indent=1)
    summary(doc["cells"])


def summary(res):
    for g, ks in GROUPS:
        rows = {k: v for k, v in res.items() if k.split("|")[0] in ks}
        if not rows:
            continue
        print(f"\n{g}: {len(rows)} cells")
        for p in PROCS:
            R = [v[p] for v in rows.values() if p in v]
            if not R:
                continue
            pen = [r["pen"] for r in R if np.isfinite(r["pen"])]
            Ks = [r["K"] for r in R if np.isfinite(r["K"])]
            print(f"  {p:<9} reached {sum(r['reached'] for r in R):>6} missed {sum(r['missed'] for r in R):>5} "
                  f"declined {sum(r['declined'] for r in R):>5} penalty {np.median(pen) if pen else np.nan:6.2f} "
                  f"runs {np.median([r['runs'] for r in R]):5.0f} K {np.median(Ks) if Ks else np.nan:7.1f}")


if __name__ == "__main__":
    main()
