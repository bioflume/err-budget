"""Declines by reason and by group (S-3.2): the decline code that run_protocol reports for
every declined search of every graded development test case, released program.

    python exp73_decline_reasons.py
"""
import collections
import os

import numpy as np

import cases as C
import pooled as P
from errbudget import run_protocol

cfg = P.CONFIGS["v1.2.1"]
total = collections.Counter()
CASES = C.all_cases()
if any(c.key == "2D-2h" for c in CASES):
    drop = "2D-2" if os.environ.get("ERRBUDGET_TABLE3_2D2", "refined") == "refined" else "2D-2h"
    CASES = [c for c in CASES if c.key != drop]
for case in CASES:
    T = case.tensor
    for ob, tgs in case.targets.items():
        for tg in tgs:
            if not case.admitted(ob, tg):
                continue
            eps = tg * abs(case.exact[ob])

            def look(i, ob=ob):
                i = tuple(i); v = T[i]
                if not np.isfinite(v[1]):
                    return np.nan, np.inf
                if not case.converged(i):
                    return np.nan, v[1]
                return v[0][ob], v[1]
            kw = dict(cfg["kw"])
            if case.floor is not None:
                kw["noise_floor"] = case.floor
            rs = collections.Counter()
            for a in case.grid():
                if not case.usable(a):
                    continue
                r = run_protocol(look, case.LAD, a, eps, hs=case.HS, kinds=case.KINDS,
                                 models=case.models, **kw)
                if not r["certified"]:
                    rs[r["declined"]] += 1
            total.update(rs)
            if rs:
                print(f"{case.key:<7} {ob:<18} {tg:<8.0e} " + ", ".join(f"{k} {v}" for k, v in sorted(rs.items())),
                      flush=True)
print("TOTAL", dict(total))
