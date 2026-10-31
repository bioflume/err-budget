"""Analysis of J3: Quantum ESPRESSO, current release against the article's pw.x 6.7.
    python analyse_j3.py [config ...]          # default: v1.2.1 v1.1
Revision 4: loads the article's QE case directly (cases.all_cases() needs SymPy for the
manufactured problems), and reports distinct accepted settings.
"""
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import qe_hpc as J
import cases as C
import pooled as P

RES = os.path.join(os.environ.get("WORK", "."), "qe_hpc")
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
R = {}
for t in J.TASKS():
    f = os.path.join(RES, t["id"], "result.json")
    if os.path.exists(f):
        R[t["id"]] = json.load(open(f))
need = ["ref", "e75", "k8", "k10"]
if not all(k in R for k in need):
    raise SystemExit("reference runs not finished")
REF = R["ref"]["obs"]
if not all(np.isfinite(R[k]["obs"][ob]) for k in need for ob in REF):
    raise SystemExit("a reference run (ref, e75, k8, k10) has no result, so nothing can be graded: "
                     + "; ".join(f"{k}: {R[k]['info']}" for k in need))
UNC = {ob: (abs(R["e75"]["obs"][ob] - REF[ob]) + abs(R["k10"]["obs"][ob] - R["k8"]["obs"][ob]))
       / abs(REF[ob]) for ob in REF}
tensor = {tuple(r["idx"]): (r["obs"], r["cost"], r["info"]) for k, r in R.items() if k.startswith("g_")}
print(f"current pw.x: {len(tensor)}/480 settings; reference {REF}; uncertainty {UNC}")
pickle.dump((tensor, REF, UNC), open(os.path.join(OUT, "qe_tensor_current.pkl"), "wb"))
old = C._qe()
oldt, OREF, _ = pickle.load(open(os.path.join(C.D, "qe_tensor.pkl"), "rb"))
common = [i for i in tensor if np.isfinite(tensor[i][1]) and np.isfinite(oldt[i][1])]
for ob in REF:
    d = [abs(tensor[i][0][ob] - oldt[i][0][ob]) / abs(OREF[ob]) for i in common]
    print(f"  {ob}: largest relative difference between releases over {len(common)} settings: "
          f"{max(d, default=float('nan')):.1e}; references differ by {abs(REF[ob]-OREF[ob])/abs(OREF[ob]):.1e}")
if len(tensor) == 480:
    new = C.Case("QEc", "Quantum ESPRESSO (current)", {i: (v[0], v[1], dict(v[2])) for i, v in tensor.items()},
                 old.LAD, old.HS, old.KINDS, dict(REF), old.targets, models=old.models, unc=UNC)
    rows = []
    for cfg in (sys.argv[1:] or ["v1.2.1", "v1.1"]):
        for case in (old, new):
            for ob, tgs in case.targets.items():
                for tg in tgs:
                    if not case.admitted(ob, tg):
                        continue
                    s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
                    s.update(config=cfg, release="6.7" if case is old else "current")
                    rows.append(s)
                    print(f"  {cfg:<6} {s['release']:<8} {ob:<6} {tg:.0e}: accepted {s['accepted']}/{s['n']} "
                          f"({s['distinct']} distinct) "
                          f"missed {s['missed']} declined(reachable) {s['declined_reach']} "
                          f"penalty {s['pen_med']:.2f}", flush=True)
    json.dump(rows, open(os.path.join(OUT, "j3_replay.json"), "w"), indent=1)
