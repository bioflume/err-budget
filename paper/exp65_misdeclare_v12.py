"""Misdeclared classes with version 1.2 (supplement, Table S2; referee point 2).

Each L3 parameter is declared as each of the two classes it does not belong to,
and the program is replayed from all 432 starts for both observables at 1e-2 and
3e-3, in the full mode (the default of version 1.2) and with the preliminary sweep.
The same is done for the four OpenFOAM 2D-1 parameters (c_D at 3e-3, c_L at 3e-2).
    python exp65_misdeclare_v12.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C
import pooled as P

CLASSES = ["discretisation", "nested-solver", "cost-only"]
NAMES = {"L3": ["N", "time steps", "Newton tol.", "GMRES tol.", "ILU drop tol."],
         "2D-1": ["n_arc", "residualControl", "GAMG tol.", "coarsest level"]}
CELLS = {"L3": [(o, t) for o in ("energy", "grad_energy") for t in (1e-2, 3e-3)],
         "2D-1": [("cd", 3e-3), ("cl", 3e-2)]}

L3 = [c for c in C._manufactured() if c.key == "L3"][0]
for case in (L3, C._dfg1()):
    cells = [(o, t) for o, t in CELLS[case.key] if o in case.targets]
    if case.key == "L3":
        cells = [(o, t) for o in case.targets for t in case.targets[o]]
    conf = [("correct", list(case.KINDS))]
    for k, nm in enumerate(NAMES[case.key]):
        for c in CLASSES:
            if c != case.KINDS[k]:
                kk = list(case.KINDS); kk[k] = c
                conf.append((f"{nm} as {c}", kk))
    for mode, scr in (("full (default)", False), ("preliminary sweep", True)):
        print("=" * 100)
        print(f"{case.label}, {mode}: cells {cells}")
        print(f"  {'declaration':<34} {'starts':>7} {'accepted':>9} {'missed':>7} "
              f"{'med pen':>8} {'worst pen':>9} {'med runs':>9}")
        for lab, kinds in conf:
            cfg = dict(P.CONFIGS["v1.2.1"]); cfg["kw"] = dict(cfg["kw"], use_screen=scr)
            cfg["kinds"] = kinds
            n = acc = miss = 0; pen = []; runs = []
            for ob, tg in cells:
                o = P.replay(case, ob, tg, cfg)
                n += o["n"]; acc += o["accepted"]; miss += o["missed"]
                pen += o["pen"]; runs += o["runs"]
            print(f"  {lab:<34} {n:>7} {acc:>9} {miss:>7} "
                  f"{np.median(pen) if pen else np.nan:>8.2f} {np.max(pen) if pen else np.nan:>9.2f} "
                  f"{np.median(runs):>9.0f}", flush=True)
