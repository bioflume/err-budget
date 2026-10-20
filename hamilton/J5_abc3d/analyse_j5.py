"""Analysis of J5, the pre-registered held-out 3D test (see PREREGISTRATION.md).
    python analyse_j5.py
The configuration and the grading rules are those fixed in PREREGISTRATION.md.
"""
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import abc_hpc as J
import cases as C
import pooled as P

RES = os.path.join(os.environ.get("WORK", "."), "abc_hpc")
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
TARGETS = {"ek": [1e-3, 1e-4], "eps": [3e-2, 1e-2, 3e-3], "umax": [1e-3, 3e-4]}

tensor = {}
for t in J.TASKS():
    f = os.path.join(RES, t["id"], "result.json")
    if not os.path.exists(f):
        continue
    r = json.load(open(f))
    if "obs" in r and r.get("cost"):
        tensor[tuple(t["idx"])] = (r["obs"], float(r["cost"]), dict(converged=True))
    else:
        tensor[tuple(t["idx"])] = ({k: np.nan for k in J.EXACT}, np.inf, dict(converged=False))
print(f"J5: {len(tensor)}/{len(J.TASKS())} settings finished, "
      f"{sum(np.isfinite(v[1]) for v in tensor.values())} usable")
pickle.dump(tensor, open(os.path.join(OUT, "j5_tensor.pkl"), "wb"))
fin = (len(J.LAD[1]) - 1, len(J.LAD[2]) - 1, len(J.LAD[3]) - 1)
print("relative error against N at the finest other settings:")
for ob in J.EXACT:
    e = [abs(tensor[(a,) + fin][0][ob] - J.EXACT[ob]) / J.EXACT[ob] if (a,) + fin in tensor
         else np.nan for a in range(len(J.LAD[0]))]
    print(f"  {ob:<5} " + "  ".join(f"{x:.2e}" for x in e))
if len(tensor) < len(J.TASKS()):
    raise SystemExit("not all runs finished")
case = C.Case("ABC3D", "OpenFOAM, 3D ABC flow (held out)", tensor, J.LAD, J.HS, J.KINDS,
              dict(J.EXACT), TARGETS, cost_note="core-seconds")
# The registered grading of the frozen program is done by paper/exp62_heldout.py, which
# reads results/j5_tensor.pkl and replays the frozen source (and, with --current, the
# released version 1.2.1).  The replay below uses the installed program, as a check.
import errbudget
print(f"replay with the installed program, errbudget {errbudget.__version__} (the registered grading is "
      "paper/exp62_heldout.py)")
rows = []
for cfg in ["v1.2.1"]:
    for ob, tgs in TARGETS.items():
        for tg in tgs:
            s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
            s["config"] = cfg
            rows.append(s)
            print(f"  {cfg:<5} {ob:<5} {tg:<7.0e} starts {s['n']} accepted {s['accepted']} missed "
                  f"{s['missed']} declined(reachable) {s['declined_reach']} declined(unreachable) "
                  f"{s['declined_unreach']} penalty {s['pen_med']:.2f}/{s['pen_max']:.2f} "
                  f"worst {s['worst']:.2f}", flush=True)
json.dump(rows, open(os.path.join(OUT, "j5_replay.json"), "w"), indent=1)
