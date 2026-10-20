"""The method on GROMACS, against GROMACS's own defaults."""
import itertools
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget import LocalBudget, run_protocol
from gromacs_water import GromacsWater

NAMES = ["rcoulomb", "fourierspacing", "pme_order", "ewald_rtol"]
LAD = [[0.9, 1.0, 1.1, 1.2, 1.4],
       [0.20, 0.16, 0.12, 0.10, 0.08],
       [4, 5, 6],
       [1e-4, 1e-5, 1e-6]]
KINDS = ["discretisation", "discretisation", "cost-only", "discretisation"]
HS = [[1.0 / c for c in LAD[0]], list(LAD[1]), [1.0 / o for o in LAD[2]],
      list(LAD[3])]
ALL = list(itertools.product(*[range(len(l)) for l in LAD]))

BINARY = os.environ.get("ERRBUDGET_GMX", "gmx_d")
CACHE = os.path.join(os.path.dirname(__file__), "data",
                     f"gromacs_tensor_{BINARY}.pkl")
if os.path.exists(CACHE):
    tensor, QS = pickle.load(open(CACHE, "rb"))
    class _P:  # the binary is needed only to build the tensor
        nmol = 1728
    p = _P()
else:
    p = GromacsWater(nbox=2, nsteps=20, binary=BINARY)
    ref, _, _, _ = p.run(p.reference())
    QS = ref["potential"]
    print(f"converged reference (beyond the ladder in every parameter): "
          f"{QS:.6f} kJ/mol")
    tensor, t0 = {}, time.time()
    for i, idx in enumerate(ALL):
        th = {NAMES[j]: LAD[j][v] for j, v in enumerate(idx)}
        obs, cost, _, info = p.run(th)
        tensor[idx] = (obs, cost)
        if i % 30 == 0:
            print(f"  {i}/{len(ALL)}  {time.time()-t0:.0f}s", flush=True)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    pickle.dump((tensor, QS), open(CACHE, "wb"))
    print(f"tensor: {len(tensor)} GROMACS runs in {time.time()-t0:.0f}s")

ok = [i for i in ALL if np.isfinite(tensor[i][0]["potential"])]
rel = lambda i: abs(tensor[i][0]["potential"] - QS) / abs(QS)
cost = lambda i: tensor[i][1]
look = lambda i: (tensor[tuple(i)][0]["potential"], tensor[tuple(i)][1])
show = lambda i: "rcoulomb=%.1f, fourier=%.2f, order=%d, rtol=%.0e" % tuple(
    LAD[k][j] for k, j in enumerate(i))

print(f"\n{'='*94}")
print(f"GROMACS [{BINARY}] {p.nmol} SPC waters ({p.nmol*3} atoms); "
      f"{len(ok)}/{len(ALL)} configurations evaluable")
print("=" * 94)

for target in [1e-5, 1e-6]:
    eps = target * abs(QS)
    feas = [i for i in ok if rel(i) <= target]
    if not feas:
        print(f"\ntarget {target:.0e}: unattainable")
        continue
    best = min(feas, key=cost)
    default = (LAD[0].index(1.0), LAD[1].index(0.12), LAD[2].index(4),
               LAD[3].index(1e-5))
    cautious = (LAD[0].index(1.4), LAD[1].index(0.08), LAD[2].index(6),
                LAD[3].index(1e-6))
    r = run_protocol(look, LAD, default, eps, hs=HS, kinds=KINDS, use_screen=False)
    rec = r["proposal"]

    print(f"\ntarget: relative error {target:.0e} in the potential energy")
    print(f"  {'':<26} {'setting':<48} {'rel. error':>11} {'ms/step':>9}")
    for lab, idx in [("GROMACS defaults", default), ("cautious user", cautious),
                     ("this method", rec), ("cheapest that attains it", best)]:
        mark = "" if rel(idx) <= target else "   MISSES"
        print(f"  {lab:<26} {show(idx):<48} {rel(idx):11.3e} "
              f"{cost(idx)*1000:8.1f}{mark}")
    print(f"  -> method vs GROMACS defaults: {cost(default)/cost(rec):.2f}x cheaper")
    print(f"  -> method vs cheapest attaining: {cost(rec)/cost(best):.2f}x, "
          f"{r['runs']} GROMACS runs used")
    if r["violations"]:
        print(f"  -> DECLARATION CONTRADICTED: "
              f"{[NAMES[k] for k in r['violations']]} declared cost-only but "
              f"moves the observable")

print(f"\n{'='*94}\nerror budget the method reports for GROMACS")
lb = LocalBudget(look, LAD, (2, 2, 1, 1), hs=HS)
lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
tot = sum(abs(x) for x in lb.rk if np.isfinite(x)) or 1.0
for k, nm in enumerate(NAMES):
    print(f"  {nm:<16} residual {lb.rk[k]:+.3e}   "
          f"{100*abs(lb.rk[k])/tot:5.1f}% of the budget")

print(f"\ncost of pme-order at rcoulomb=1.0, fourier=0.12, rtol=1e-5 "
      f"(GROMACS default order is 4):")
for j, o in enumerate(LAD[2]):
    i = (LAD[0].index(1.0), LAD[1].index(0.12), j, LAD[3].index(1e-5))
    print(f"  order {o}:  rel. error {rel(i):.3e}   {cost(i)*1000:.1f} ms/step")
