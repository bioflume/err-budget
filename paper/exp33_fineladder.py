"""A test designed to break the method.

Section 6.4 concluded that what carries the reliability is not the accuracy of
the error estimate -- which is optimistic at a fifth of anchors -- but the
coarseness of the ladder: rungs differing by a factor of two in the mesh and a
hundred in a tolerance mean that a configuration passing the test usually sits
well inside the target rather than at its boundary, so a 20% optimistic estimate
rarely changes a decision.

That account makes a prediction: on a FINELY spaced ladder the protection erodes
and reliability should fall.  Predictions of this kind are worth more tested than
stated, so here the same problem is given ladders whose consecutive rungs differ
by about 20% in the discretisation parameters and a factor of three in the
tolerances, instead of factors of two and a hundred.  Everything else --- the
method, the margin, the protocol --- is unchanged.

If safety holds up, the explanation in 6.4 is wrong.  If it falls, the paper has
identified the condition under which the method should not be trusted, which is
more useful than another passing test.
"""
import itertools
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget import run_protocol
from errbudget.problems import Burgers1D

D = os.path.join(os.path.dirname(__file__), "data")
NAMES = ["N", "steps", "newton_tol", "gmres_tol"]

# ~20% between rungs in the mesh and the step, ~3x in the tolerances, against
# 2x and 100x in the coarse ladder of the main study.
FINE = [[96, 116, 140, 168, 202],
        [70, 84, 101, 121, 145],
        [1e-4, 3e-5, 1e-5, 3e-6, 1e-6],
        [1e-4, 3e-5, 1e-5, 3e-6, 1e-6]]
COARSE = [[32, 64, 128, 256, 512],
          [25, 50, 100, 200, 400],
          [1e-2, 1e-4, 1e-6, 1e-8, 1e-10],
          [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
KIN = ["discretisation", "discretisation", "nested-solver", "nested-solver"]

p = Burgers1D(nu=0.05, T=1.0)
EX = p.exact_observables()
CACHE = os.path.join(D, "burgers_fine_tensor.pkl")
ALL = list(itertools.product(*[range(5) for _ in range(4)]))

if os.path.exists(CACHE):
    tensor = pickle.load(open(CACHE, "rb"))
else:
    tensor, t0 = {}, time.time()
    for i, idx in enumerate(ALL):
        th = {NAMES[j]: FINE[j][v] for j, v in enumerate(idx)}
        obs, C, w, _ = p.run(th)
        tensor[idx] = (obs, C)
        if i % 50 == 0:
            print(f"  {i}/{len(ALL)}  {time.time()-t0:.0f}s", flush=True)
    pickle.dump(tensor, open(CACHE, "wb"))
    print(f"fine tensor: {len(tensor)} runs in {time.time()-t0:.0f}s", flush=True)

HSF = [[2 * np.pi / n for n in FINE[0]], [1.0 / s for s in FINE[1]],
       list(FINE[2]), list(FINE[3])]


def grade(tensor, LAD, HS, ob, targets_rel):
    QS = EX[ob]
    err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
    lk = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    out = []
    for tr in targets_rel:
        eps = tr * abs(QS)
        feas = [i for i in ALL if err[i] <= eps]
        if not feas or len(feas) == len(ALL):
            out.append((tr, None, None, None, None, len(feas)))
            continue
        c_opt = min(tensor[i][1] for i in feas)
        cert = safe = 0
        pens, slack = [], []
        for a in ALL:
            r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KIN, use_screen=False)
            if not r["certified"]:
                continue
            cert += 1
            q = r["proposal"]
            safe += err[q] <= eps
            pens.append(tensor[q][1] / c_opt)
            slack.append(err[q] / eps)      # how far inside the target it lands
        out.append((tr, cert, safe, float(np.median(pens)),
                    float(np.median(slack)), len(feas)))
    return out


print("=" * 100)
print("the same method, the same problem, ladders ~20% apart instead of 2x apart")
print("=" * 100)
print(f"  coarse ladder mesh   : {COARSE[0]}   (2.0x between rungs)")
print(f"  fine ladder mesh     : {FINE[0]}   (1.2x between rungs)")
print(f"  coarse ladder tol    : 1e-2 ... 1e-10  (100x between rungs)")
print(f"  fine ladder tol      : 1e-4 ... 1e-6   (3x between rungs)")
print()
print(f"  {'observable':<14} {'target':>9} {'certified':>10} {'attained':>9} "
      f"{'safety':>8} {'med penalty':>12} {'median slack':>13}")
tot_c = tot_s = 0
for ob in ["energy", "dissipation"]:
    for tr, cert, safe, pen, slack, nfeas in grade(
            tensor, FINE, HSF, ob, [3e-3, 1e-3, 3e-4, 1e-4]):
        if cert is None:
            print(f"  {ob:<14} {tr:>9.0e} "
                  f"{'no discrimination on this grid':>46} ({nfeas}/625 attain)")
            continue
        tot_c += cert; tot_s += safe
        print(f"  {ob:<14} {tr:>9.0e} {cert:>10} {safe:>9} "
              f"{100*safe/max(cert,1):>7.1f}% {pen:>11.2f}x {slack:>12.3f}")
print()
print(f"  overall on the fine ladder: {tot_s}/{tot_c} certified recommendations "
      f"attained the target ({100*tot_s/max(tot_c,1):.2f}%)")
print()
print("  'median slack' is the true error of the recommendation as a fraction of")
print("  the target.  On the coarse ladder it is small because a rung overshoots;")
print("  the prediction of 6.4 is that it approaches 1 here, leaving no room for")
print("  an optimistic estimate to be absorbed.")
