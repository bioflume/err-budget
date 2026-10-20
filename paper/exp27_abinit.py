"""A third production code: ABINIT, plane-wave DFT on bulk silicon.

This code serves the third audience the paper addresses -- electronic structure
-- and it tests a prediction made before the run.  A Quantum ESPRESSO developer,
answering a user on the PWSCF mailing list, wrote that a convergence threshold
may be "acceptable if you just want the energy but not if you want forces".
That is the paper's observable-dependence finding, stated by a code author as
folklore.  If it is a property of plane-wave DFT rather than a quirk of one
code, the method must recover it here, from the outside, without being told.

Parameters are ABINIT's own and map onto the taxonomy exactly as the
applicability table predicts:

    ecut    plane-wave cutoff         discretisation
    ngkpt   Brillouin-zone sampling   discretisation
    toldfe  SCF convergence threshold nested solver
    nline   line minimisations/step   cost-only

Ground truth is a reference set strictly beyond the ladder in every parameter,
so no tested configuration inherits a small error by sharing a value with it.
"""
import itertools
import json
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget import LocalBudget, run_protocol
from abinit_silicon import AbinitSilicon

NAMES = ["ecut", "ngkpt", "toldfe", "nline"]
LAD = [[8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0], [2, 3, 4, 6],
       [1e-6, 1e-8, 1e-10, 1e-12], [2, 4, 6]]
KINDS = ["discretisation", "discretisation", "nested-solver", "cost-only"]
HS = [[1.0 / np.sqrt(e) for e in LAD[0]],      # resolution ~ 1/sqrt(ecut)
      [1.0 / k for k in LAD[1]], list(LAD[2]), [1.0 / n for n in LAD[3]]]
# Version 1.1: a plane-wave cutoff converges exponentially, not algebraically
# in 1/sqrt(ecut).  The exponential model fits e(x) = c (exp(-b x) - exp(-b x_a))
# in the cutoff itself.
HS_EXP = [list(LAD[0])] + HS[1:]
# The exponential model extrapolates further than the power law; its own fit
# residual sets the margin (the residual-gated margin of Section 7.4, alpha = 5).
MODELS = {"power": (HS, None, {}),
          "exponential": (HS_EXP, ["exponential", "power", "power", "power"], {}),
          "exp. + gated": (HS_EXP, ["exponential", "power", "power", "power"],
                           dict(resid_margin=5.0))}
ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
CACHE = os.path.join(os.path.dirname(__file__), "data", "abinit_tensor_v11.pkl")

p = AbinitSilicon()
if os.path.exists(CACHE):
    tensor, REF = pickle.load(open(CACHE, "rb"))
else:
    print("converged reference (beyond the ladder in every parameter) ...",
          flush=True)
    t0 = time.time()
    REF, _, _, _ = p.run(p.reference())
    print(f"  E = {REF['energy']:.10f} Ha, |F| = {REF['force']:.8e} Ha/Bohr "
          f"({time.time()-t0:.0f}s)", flush=True)
    tensor, t0 = {}, time.time()
    for i, idx in enumerate(ALL):
        th = {NAMES[j]: LAD[j][v] for j, v in enumerate(idx)}
        obs, cost, wall, info = p.run(th)
        tensor[idx] = (obs, cost)
        if i % 16 == 0:
            print(f"  {i}/{len(ALL)}  {time.time()-t0:.0f}s", flush=True)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    pickle.dump((tensor, REF), open(CACHE, "wb"))
    print(f"tensor: {len(tensor)} ABINIT runs in {time.time()-t0:.0f}s", flush=True)

ok = [i for i in ALL if np.isfinite(tensor[i][0]["energy"])]
show = lambda i: "ecut=%g, ngkpt=%d, toldfe=%.0e, nline=%d" % tuple(
    LAD[k][j] for k, j in enumerate(i))

print(f"\n{'='*100}")
print(f"ABINIT 9.10.4, bulk Si with one displaced atom; "
      f"{len(ok)}/{len(ALL)} configurations evaluable")
print(f"reference: E = {REF['energy']:.10f} Ha, |F| = {REF['force']:.8e} Ha/Bohr")
print("=" * 100)

out = {}
TARGETS = {"energy": [1e-2, 3e-3, 1e-3], "force": [1e-1, 3e-2, 1e-2]}
START = (1, 1, 1, 1)   # the ABINIT tutorial's silicon settings: ecut 12, 3x3x3, 1e-8, 4
for ob in ["energy", "force"]:
    QS = REF[ob]
    rel = lambda i, ob=ob, QS=QS: abs(tensor[i][0][ob] - QS) / abs(QS)
    look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    print(f"\n--- observable: {ob} "
          f"(reference {QS:.8g}) {'-'*(60-len(ob))}")
    for target in TARGETS[ob]:
        eps = target * abs(QS)
        feas = [i for i in ok if rel(i) <= target]
        if not feas:
            print(f"  target {target:.0e} rel: unattainable on the grid "
                  f"(best {min(rel(i) for i in ok):.2e})")
            continue
        best = min(feas, key=lambda i: tensor[i][1])
        print(f"  target {target:.0e} rel")
        print(f"    {'tutorial-like start':<26} {show(START):<50} "
              f"{rel(START):.2e}  {tensor[START][1]:6.1f}s"
              f"{'' if rel(START)<=target else '   MISSES'}")
        for mname, (hs, models, kw) in MODELS.items():
            r = run_protocol(look, LAD, START, eps, hs=hs, kinds=KINDS,
                             use_screen=False, models=models, **kw)
            q = r["proposal"]
            lab = f"method, {mname} model"
            print(f"    {lab:<26} {show(q):<50} "
                  f"{rel(q):.2e}  {tensor[q][1]:6.1f}s"
                  f"{'' if r['certified'] else '   (declined)'}"
                  f"{'   MISSES' if r['certified'] and rel(q) > target else ''}"
                  f"   est {abs(r['budget'].E0)/abs(QS):.2e}, {r['runs']} runs")
            out[(ob, target, mname)] = dict(proposal=q, certified=r["certified"],
                                            rel=rel(q), cost=tensor[q][1])
        print(f"    {'cheapest attaining':<26} {show(best):<50} "
              f"{rel(best):.2e}  {tensor[best][1]:6.1f}s")
        # every start, as in the exhaustive replays of Section 7
        for mname, (hs, models, kw) in MODELS.items():
            cert = miss = 0; pen, runs = [], []
            for a in ok:
                r = run_protocol(look, LAD, a, eps, hs=hs, kinds=KINDS,
                                 use_screen=False, models=models, **kw)
                runs.append(r["runs"])
                if r["certified"]:
                    cert += 1; q = r["proposal"]
                    miss += rel(q) > target
                    pen.append(tensor[q][1] / tensor[best][1])
            print(f"    all {len(ok)} starts, {mname:<11} model: {cert} certified, "
                  f"{miss} missed, median penalty "
                  f"{np.median(pen) if pen else float('nan'):.2f}, median runs "
                  f"{np.median(runs):.0f}")

# ----------------------------------------------------- the pre-registered test
print(f"\n{'='*100}")
PRED = 1e-2
print("the prediction: a setting acceptable for the energy is not acceptable")
print("for the force.  Both observables come from the SAME runs.")
print("=" * 100)
for ob in ["energy", "force"]:
    QS = REF[ob]
    rel = lambda i, QS=QS, ob=ob: abs(tensor[i][0][ob] - QS) / abs(QS)
    feas = [i for i in ok if rel(i) <= PRED]
    if feas:
        b = min(feas, key=lambda i: tensor[i][1])
        print(f"  cheapest setting reaching {PRED:.0e} relative in the {ob:<7}: "
              f"{show(b)}   ({tensor[b][1]:.1f}s)")
for ob, other in [("energy", "force"), ("force", "energy")]:
    QS, QO = REF[ob], REF[other]
    rel = lambda i, QS=QS, ob=ob: abs(tensor[i][0][ob] - QS) / abs(QS)
    relo = lambda i, QO=QO, other=other: abs(tensor[i][0][other] - QO) / abs(QO)
    feas = [i for i in ok if rel(i) <= PRED]
    if not feas:
        continue
    b = min(feas, key=lambda i: tensor[i][1])
    print(f"  that setting gives {relo(b):.2e} relative in the {other} "
          f"-- {'PASSES' if relo(b) <= PRED else 'FAILS'} the same {PRED:.0e} target")

print(f"\n{'='*100}\nerror budget the method reports for ABINIT")
for ob in ["energy", "force"]:
    look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    for mname, (hs, models, kw) in MODELS.items():
      lb = LocalBudget(look, LAD, START, hs=hs, models=models, **kw)
      lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
      tot = sum(abs(x) for x in lb.rk if np.isfinite(x)) or 1.0
      true = tensor[START][0][ob] - REF[ob]
      print(f"  {ob}, {mname} model at the tutorial start: estimate "
            f"{sum(x for x in lb.rk if np.isfinite(x)):+.3e}, true error {true:+.3e}")
      for k, nm in enumerate(NAMES):
        print(f"    {nm:<8} residual {lb.rk[k]:+.3e}   "
              f"{100*abs(lb.rk[k])/tot:5.1f}% of the budget   "
              f"(declared {KINDS[k]})")
