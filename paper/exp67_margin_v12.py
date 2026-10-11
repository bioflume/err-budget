"""Choice of the base safety factor of version 1.2 (Appendix B).

The base factor F (fs_cap; the gated factor is F (1 + 5 q)) is scanned over every
admitted cell of the development problems with exact solutions, L2-L5 and 2D-1,
with the rest of version 1.2 unchanged.  For each F: accepted, missed and the median
of the per-cell median penalties, per problem; and the smallest F that is safe when
each problem in turn is the calibration set.
    python exp67_margin_v12.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C
import pooled as P

FS = [1.0, 1.1, 1.2, 1.25, 1.35, 1.5, 2.0]
cs = [c for c in C._manufactured()] + [C._dfg1()]
res = {}
for F in FS:
    cfg = dict(P.CONFIGS["v1.2.1"]); cfg["kw"] = dict(cfg["kw"], fs_cap=F)
    for case in cs:
        acc = miss = 0; pm = []
        for ob, tgs in case.targets.items():
            for tg in tgs:
                if not case.admitted(ob, tg):
                    continue
                s = P.summarise(P.replay(case, ob, tg, cfg))
                acc += s["accepted"]; miss += s["missed"]
                if np.isfinite(s["pen_med"]):
                    pm.append(s["pen_med"])
        res[(F, case.key)] = (acc, miss, float(np.median(pm)) if pm else np.nan)
        print(f"F={F:<5} {case.key:<5} accepted {acc:<6} missed {miss:<4} median penalty "
              f"{res[(F, case.key)][2]:.2f}", flush=True)
print("\nsmallest safe F when calibrated on one problem, and the misses it leaves elsewhere:")
for key in [c.key for c in cs]:
    safe = [F for F in FS if res[(F, key)][1] == 0]
    F0 = min(safe) if safe else None
    if F0 is None:
        print(f"  {key}: no F in the scan is safe"); continue
    other = {k2: res[(F0, k2)][1] for k2 in [c.key for c in cs] if k2 != key}
    print(f"  calibrated on {key}: F = {F0}; misses elsewhere {other}")
