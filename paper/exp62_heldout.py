"""Held-out tests H1 and H2 (and H3 if its Hamilton results are present), replayed
with the frozen program exactly as fixed in PREREGISTRATION.md.

The script refuses to run if the program's source differs from the frozen one.
    python exp62_heldout.py
"""
import glob
import hashlib
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cases as C
import pooled as P

FROZEN = os.environ.get("KNOBS_FROZEN_HASH", "")
SRC = os.path.join(HERE, "..", "src", "knobs")


def source_hash():
    sums = "".join(hashlib.sha256(open(f, "rb").read()).hexdigest() + "  " + os.path.basename(f) + "\n"
                   for f in sorted(glob.glob(os.path.join(SRC, "*.py"))))
    return hashlib.sha256(sums.encode()).hexdigest()


reg = open(os.path.join(HERE, "..", "PREREGISTRATION.md")).read()
h = source_hash()
if h not in reg:
    raise SystemExit(f"source hash {h} is not the frozen one recorded in PREREGISTRATION.md")
print(f"frozen program verified (sha256 {h[:16]}...)")


def heldout_case(key):
    if key == "H1":
        import heldout_fem as M
        exact, tg = M.EXACT, {"mean": [1e-3, 3e-4], "flux": [3e-2, 1e-2, 3e-3]}
        label = "H1: 3D Poisson, scikit-fem"
    else:
        import heldout_cavity as M
        exact, tg = M.REF, {o: [3e-2, 1e-2, 3e-3] for o in M.REF}
        label = "H2: cavity Re 1000, OpenFOAM"
    t = pickle.load(open(os.path.join(C.D, f"heldout_{key}.pkl"), "rb"))
    t = {i: (v[0], v[1], dict(v[2])) for i, v in t.items()}
    return C.Case(key, label, t, M.LAD, M.HS, M.KINDS, dict(exact), tg)


rows = []
cases = [heldout_case("H1"), heldout_case("H2")]
j5 = os.path.join(HERE, "..", "hamilton", "results", "j5_tensor.pkl")
if os.path.exists(j5):
    sys.path.insert(0, os.path.join(HERE, "..", "hamilton", "J5_abc3d"))
    sys.path.insert(0, os.path.join(HERE, "..", "hamilton", "common"))
    import abc_hpc as A
    t = pickle.load(open(j5, "rb"))
    cases.append(C.Case("H3", "H3: 3D ABC flow, OpenFOAM (Hamilton)", t, A.LAD, A.HS, A.KINDS,
                        dict(A.EXACT), {"ek": [1e-3, 1e-4], "eps": [3e-2, 1e-2, 3e-3],
                                        "umax": [1e-3, 3e-4]}))
for case in cases:
    G = case.grid()
    print(f"\n{case.label}: {len(G)} settings, {sum(case.usable(i) for i in G)} usable, "
          f"{sum(case.usable(i) and not case.converged(i) for i in G)} not converged")
    fin = [len(l) - 1 for l in case.LAD]
    for ob in case.targets:
        e = []
        for a in range(len(case.LAD[0])):
            i = tuple([a] + fin[1:])
            e.append(abs(case.tensor[i][0][ob] - case.exact[ob]) / abs(case.exact[ob])
                     if case.usable(i) else np.nan)
        print(f"  {ob:<6} error against the first parameter at the finest others: "
              + "  ".join(f"{x:.2e}" for x in e))
    for cfg in ("v1.2", "v1.1"):
        for ob, tgs in case.targets.items():
            for tg in tgs:
                s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
                s.update(config=cfg, case=case.key)
                rows.append(s)
                print(f"  {cfg:<5} {ob:<6} {tg:<7.0e} starts {s['n']:<4} start-failed "
                      f"{s['start_failed']:<4} accepted {s['accepted']:<4} missed {s['missed']:<3} "
                      f"declined(reachable) {s['declined_reach']:<4} declined(unreachable) "
                      f"{s['declined_unreach']:<4} penalty {s['pen_med']:.2f}/{s['pen_max']:.2f} "
                      f"runs {s['runs_med']:.0f} worst {s['worst']:.2f}", flush=True)
json.dump(rows, open(os.path.join(C.D, "heldout_replay.json"), "w"), indent=1)
