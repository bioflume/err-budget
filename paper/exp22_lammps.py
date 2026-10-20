"""Running the method on LAMMPS, and comparing with LAMMPS's own guidance.

LAMMPS documents `kspace_style pppm <accuracy>` as setting "the relative RMS
error in per-atom forces", used "in conjunction with the pairwise cutoff to
determine ... the grid size".  A user following that documentation sets the
accuracy to the precision they want and chooses the cutoff by convention; the
interpolation order is left at its default of 5.  We compare that practice with
what the method recommends when it is given the same three parameters and the
same accuracy target.

Ground truth is the Madelung constant, so the error is known exactly.
"""
import itertools
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget import run_protocol
from lammps_madelung import LammpsMadelung
from lammps_subproc import run_one

NAMES = ["cutoff", "kspace_accuracy", "kspace_order"]
LAD = ([[6.0, 8.0, 10.0, 12.0, 14.0], [1e-2, 1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]]
       if int(os.environ.get("ERRBUDGET_LMP_N", "20")) <= 20 else
       [[8.0, 10.0, 12.0, 14.0], [1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]])
KINDS = ["discretisation", "discretisation", "cost-only"]
# refinement scalars: cutoff refines by growing, accuracy by shrinking, and the
# interpolation order is accuracy-neutral so its scalar is nominal
HS = [[1.0 / c for c in LAD[0]], list(LAD[1]), [1.0 / o for o in LAD[2]]]
ALL = list(itertools.product(*[range(len(l)) for l in LAD]))

NCELLS = int(os.environ.get("ERRBUDGET_LMP_N", "20"))
p = LammpsMadelung(n_cells=NCELLS, n_eval=3)
EXACT = p.exact_observables()
OB = "energy_per_ion"
QS = EXACT[OB]

CACHE = os.path.join(os.path.dirname(__file__), "data",
                     f"lammps_tensor_n{NCELLS}.pkl")
if os.path.exists(CACHE):
    tensor = pickle.load(open(CACHE, "rb"))
else:
    tensor, t0 = {}, time.time()
    nbad = 0
    for i, idx in enumerate(ALL):
        th = {NAMES[j]: LAD[j][v] for j, v in enumerate(idx)}
        obs, cost, info = run_one(NCELLS, th, n_eval=3)
        tensor[idx] = (obs, cost)
        if info.get("unavailable"):
            nbad += 1
            print(f"  {th} UNAVAILABLE: {info['unavailable']}", flush=True)
        if i % 8 == 0:
            print(f"  {i}/{len(ALL)}  {time.time()-t0:.0f}s  "
                  f"({nbad} unavailable)", flush=True)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    pickle.dump(tensor, open(CACHE, "wb"))
    print(f"tensor: {len(tensor)} LAMMPS runs in {time.time()-t0:.0f}s")

look = lambda i: (tensor[tuple(i)][0][OB], tensor[tuple(i)][1])
rel = lambda i: abs(tensor[i][0][OB] - QS) / abs(QS)
cost = lambda i: tensor[i][1]
show = lambda i: "cutoff=%.0f, accuracy=%.0e, order=%d" % tuple(
    LAD[k][j] for k, j in enumerate(i))

print(f"\n{'='*92}")
print(f"LAMMPS {p.natoms}-ion rocksalt, Madelung ground truth "
      f"({QS:.6f} eV/ion), {len(ALL)} configurations")
print("=" * 92)

for target in [1e-4, 1e-5]:
    eps = target * abs(QS)
    feas = [i for i in ALL if rel(i) <= target]
    if not feas:
        print(f"\ntarget {target:.0e}: unattainable within these values")
        continue
    best = min(feas, key=cost)

    # (A) what LAMMPS's documentation leads a user to do
    a = (LAD[0].index(10.0), LAD[1].index(target), LAD[2].index(5))
    # (B) a cautious user, tightening everything
    tight = min(LAD[1], key=lambda v: abs(np.log10(v) - np.log10(target * 1e-2)))
    b = (LAD[0].index(14.0), LAD[1].index(tight), LAD[2].index(7))
    # (C) the method, started from the documented practice
    r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, use_screen=False)
    c = r["proposal"]

    print(f"\ntarget: relative error {target:.0e} in the energy per ion")
    print(f"  {'':<26} {'setting':<44} {'rel. error':>11} {'cost/eval':>10}")
    for lab, idx in [("LAMMPS documented practice", a), ("cautious user", b),
                     ("this method", c), ("cheapest that attains it", best)]:
        mark = "" if rel(idx) <= target else "   MISSES TARGET"
        print(f"  {lab:<26} {show(idx):<44} {rel(idx):11.3e} {cost(idx):9.4f}s{mark}")
    print(f"  -> method vs documented practice: {cost(a)/cost(c):.2f}x cheaper, "
          f"error {rel(c)/rel(a):.2f}x")
    print(f"  -> method vs cheapest attaining:  {cost(c)/cost(best):.2f}x, "
          f"using {r['runs']} LAMMPS runs to find it")

# ---------------------------------------------------------------- diagnostics
print(f"\n{'='*92}\nwhat the method reports about each LAMMPS parameter")
from errbudget import LocalBudget
lb = LocalBudget(look, LAD, (2, 2, 2), hs=HS)
lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
tot = sum(abs(x) for x in lb.rk if np.isfinite(x)) or 1.0
for k, nm in enumerate(NAMES):
    print(f"  {nm:<18} residual {lb.rk[k]:+.3e}   {100*abs(lb.rk[k])/tot:5.1f}% of the budget")

print(f"\ncost of the interpolation order at cutoff=10, accuracy=1e-5 "
      f"(LAMMPS adjusts the grid to hold accuracy fixed):")
for j, o in enumerate(LAD[2]):
    i = (LAD[0].index(10.0), LAD[1].index(1e-5), j)
    print(f"  order {o}:  rel. error {rel(i):.3e}   cost {cost(i):.4f}s")
