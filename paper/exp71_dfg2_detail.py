"""DFG 2D-2: every target, including those below ten times the reference
uncertainty (not admitted to the pooled table), for the configurations that
differ in how interacting errors are treated (referee points 4 and 8).
    python exp71_dfg2_detail.py
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C
import pooled as P

case = C._dfg2()
print("reference uncertainty (relative):", {k: f"{v:.1e}" for k, v in case.unc.items()})
for ob, tgs in case.targets.items():
    for tg in tgs:
        adm = "admitted" if case.admitted(ob, tg) else "NOT admitted"
        for cfg in ("v1.1", "v1.1+gated", "v1.2-credit", "v1.2-nojoint", "v1.2", "v1.2.1"):
            s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
            print(f"{ob:<7} {tg:<7.0e} {adm:<13} {cfg:<13} starts {s['n']} start-failed {s['start_failed']} "
                  f"accepted {s['accepted']} missed {s['missed']} declined(reach) {s['declined_reach']} "
                  f"pen {s['pen_med']:.2f} worst {s['worst']:.2f}", flush=True)
