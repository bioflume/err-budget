"""Item 7: reliability with eight parameters, half of them nearly inert.

L3 with the three further solver parameters of the class-inference study
(Section 3.1): gmres_restart, precond_freq (how many time steps a preconditioner
is reused) and newton_max.  The 8-parameter grid has 11,664 settings and is
not built; the method is run live from 200 random starts at one target
(energy, 1e-2 relative) and every recommendation is graded against the exact
solution.  Runs are cached and shared between starts, and the calibration cost
of each start counts only the runs that start requested.
"""
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from knobs import run_protocol
from knobs.problems import Burgers2D
import grids

NAMES = ["N", "steps", "newton_tol", "gmres_tol", "ilu_drop", "gmres_restart",
         "precond_freq", "newton_max"]
LAD = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8], [1e-2, 1e-5, 1e-8],
       [1e-1, 1e-2, 1e-3], [10, 30, 100], [20, 5, 1], [2, 5, 20]]
KINDS = ["discretisation", "discretisation", "nested-solver", "nested-solver",
         "cost-only", "nested-solver", "cost-only", "nested-solver"]
HS = [[2 * np.pi / n for n in LAD[0]], [1.0 / s for s in LAD[1]], list(LAD[2]),
      list(LAD[3]), list(LAD[4]), [1.0 / r for r in LAD[5]], list(LAD[6]),
      [1.0 / m for m in LAD[7]]]
OB, TG = "energy", 1e-2
NSTART = int(os.environ.get("KNOBS_NSTART", 200))
CACHE = os.path.join(os.path.dirname(__file__), "data", "eight_cache.pkl")

b = Burgers2D(nu=0.05, T=1.0)
EX = b.exact_observables()
eps = TG * abs(EX[OB])
store = pickle.load(open(CACHE, "rb")) if os.path.exists(CACHE) else {}
t0, new = time.time(), [0]


def look(idx):
    idx = tuple(idx)
    if idx not in store:
        th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
        obs, c, w, info = b.run(th)
        store[idx] = (obs, c)
        new[0] += 1
        if new[0] % 50 == 0:
            pickle.dump(store, open(CACHE, "wb"))
            print(f"    {len(store)} distinct runs, {time.time()-t0:.0f}s", flush=True)
    return store[idx][0][OB], store[idx][1]


rng = np.random.default_rng(2026)
starts = sorted({tuple(int(rng.integers(0, len(l))) for l in LAD) for _ in range(3 * NSTART)})
starts = [starts[i] for i in rng.permutation(len(starts))[:NSTART]]
P5 = grids.L3
t5 = grids.load(P5)
feas5 = [i for i in grids.grid(P5) if abs(t5[i][0][OB] - EX[OB]) <= eps]
c5 = min(t5[i][1] for i in feas5)

for mode, scr in [("screened (default)", True), ("full", False)]:
    cert = miss = 0; pens, runs, calib, recs = [], [], [], set(); worst = 0.0
    for n, a in enumerate(starts):
        r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, use_screen=scr)
        runs.append(r["runs"])
        calib.append(r["protocol_cost"] / store[a][1])
        if r["certified"]:
            q = r["proposal"]; cert += 1; recs.add(q)
            e = abs(store[q][0][OB] - EX[OB]) / eps
            worst = max(worst, e); miss += e > 1
            pens.append(store[q][1] / c5)
        if n % 20 == 19:
            print(f"  {mode}: {n+1}/{len(starts)} starts done", flush=True)
    pickle.dump(store, open(CACHE, "wb"))
    seen = [i for i in store if abs(store[i][0][OB] - EX[OB]) <= eps]
    cbest = min(store[i][1] for i in seen)
    print(f"\n8-parameter L3, {mode} mode, energy 1e-2, {len(starts)} starts:")
    print(f"  certified {cert}, missed {miss}, worst error {worst:.2f} eps, "
          f"distinct recommendations {len(recs)}")
    print(f"  median runs {np.median(runs):.0f} (max {max(runs)}), median calibration cost "
          f"{np.median(calib):.0f} production runs at the start setting")
    print(f"  penalty against the 5-parameter grid optimum: median {np.median(pens):.2f}, "
          f"worst {max(pens):.2f}; against the cheapest attaining setting among the "
          f"{len(store)} 8-parameter runs made: median "
          f"{np.median([p * c5 / cbest for p in pens]):.2f}")
