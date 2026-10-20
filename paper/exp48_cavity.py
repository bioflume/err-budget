"""L5: driven cavity with a boundary layer; non-integral observables.

Builds the exhaustive tensor (7 x 5 x 5 x 4 = 700 runs; resumable), then
replays the method from every start, with and without the residual-gated
margin, and the baselines, and reports the observed convergence orders.
"""
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget.problems import DrivenCavity
from replay import grid, method_cell, baseline_cell, fmt_method
import grids

NAMES, LAD, KINDS = grids.L5["names"], grids.L5["lad"], grids.L5["kinds"]
HS = grids.hs(grids.L5)
TARGETS = grids.L5["targets"]
CACHE = os.path.join(os.path.dirname(__file__), "data", "cavity_tensor.pkl")

p = DrivenCavity()
EX = p.exact_observables()
if os.path.exists(CACHE):
    tensor = pickle.load(open(CACHE, "rb"))
else:
    part = CACHE + ".partial"
    tensor = pickle.load(open(part, "rb")) if os.path.exists(part) else {}
    t0 = time.time()
    todo = [i for i in grid(LAD) if i not in tensor]
    for n, idx in enumerate(todo):
        th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
        obs, C, w, info = p.run(th)
        tensor[idx] = (obs, C, info)
        if n % 20 == 0 or n == len(todo) - 1:
            pickle.dump(tensor, open(part, "wb"))
            print(f"  {len(tensor)}/{len(grid(LAD))}  {time.time()-t0:.0f}s", flush=True)
    pickle.dump(tensor, open(CACHE, "wb"))
    os.remove(part)
    print(f"tensor: {len(tensor)} runs in {time.time()-t0:.0f}s")

print("=" * 100)
print("L5  driven cavity with a boundary layer (delta = 0.05, nu = 0.01), "
      f"{len(tensor)} configurations")
print("exact values:", {k: round(v, 10) for k, v in EX.items()})
print("=" * 100)

# observed order of convergence at the tightest tolerances
fin = (len(LAD[1]) - 1, len(LAD[2]) - 1, len(LAD[3]) - 1)
print("\nerror against N at the tightest tolerances (observed order between levels):")
for ob in EX:
    e = [abs(tensor[(j,) + fin][0][ob] - EX[ob]) / abs(EX[ob]) for j in range(len(LAD[0]))]
    ords = [np.log(e[j] / e[j + 1]) / np.log(LAD[0][j + 1] / LAD[0][j]) for j in range(len(e) - 1)]
    print(f"  {ob:<10} " + "  ".join(f"{x:.1e}" for x in e))
    print(f"  {'order':<10} " + "  ".join(f"{x:7.2f}" for x in ords))

print("\nshare of the configurations whose Newton iteration did not converge: "
      f"{sum(not v[2]['converged'] for v in tensor.values())}/{len(tensor)}")

for ob, tgs in TARGETS.items():
    for tg in tgs:
        eps = tg * abs(EX[ob])
        print(f"\n--- {ob}, target {tg:.0e} relative")
        for lab, kw in [("method", {}), ("method, gated margin alpha=5", dict(resid_margin=5.0))]:
            o = method_cell(tensor, LAD, HS, KINDS, ob, EX[ob], eps, **kw)
            print(f"  {lab:<30} " + fmt_method(o))
            for m in o["misses"][:5]:
                print(f"      miss: start {m[0]} -> {m[1]}, true {m[2]:.2f} eps, estimate {m[3]:.2f} eps")
        for lab, r in baseline_cell(tensor, LAD, ob, EX[ob], eps).items():
            print(f"  {lab:<30} attained {r['attained']}/{r['n']}, penalty median "
                  f"{r['pen']:.2f}, runs median {r['runs']:.0f}")
