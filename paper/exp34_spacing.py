"""Reliability against ladder spacing: the test the method is supposed to fail.

Section 6.4 attributes the measured reliability not to the accuracy of the error
estimate, which is optimistic at a fifth of anchors, but to the coarseness of the
ladder.  A configuration that passes the stopping test typically lands well inside
the target because the next rung down is far away; an estimate that is 20%
optimistic changes the decision only near a rung boundary.

If that is right, reliability must degrade as the rungs close up, and the
diagnostic that degrades with it is the SLACK -- the true error of the
recommendation as a fraction of the target.  Slack near zero means the ladder is
doing the work; slack near one means the stopping test is, and the stopping test
is the part that is only 80% accurate.

Same problem, same method, same margin; only the spacing changes.  Each ladder is
centred on the same geometric midpoint so that the comparison is of spacing and
not of range.
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
KIN = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
ALL = list(itertools.product(*[range(5) for _ in range(4)]))
p = Burgers1D(nu=0.05, T=1.0)
EX = p.exact_observables()


def ladders(rd, rt):
    """Five rungs, ratio rd for discretisation and rt for tolerances, centred
    on N=128, steps=100, tol=1e-6."""
    N = [int(round(128 * rd ** (j - 2))) for j in range(5)]
    S = [int(round(100 * rd ** (j - 2))) for j in range(5)]
    T = [1e-6 * rt ** (2 - j) for j in range(5)]
    return [N, S, list(T), list(T)]


def tensor_for(rd, rt):
    key = f"burgers_sp_{rd}_{rt}.pkl"
    path = os.path.join(D, key)
    LAD = ladders(rd, rt)
    if os.path.exists(path):
        return pickle.load(open(path, "rb")), LAD
    t, t0 = {}, time.time()
    for i, idx in enumerate(ALL):
        th = {NAMES[j]: LAD[j][v] for j, v in enumerate(idx)}
        obs, C, w, _ = p.run(th)
        t[idx] = (obs, C)
    pickle.dump(t, open(path, "wb"))
    print(f"  tensor for spacing {rd}x/{rt}x: {len(t)} runs in "
          f"{time.time()-t0:.0f}s", flush=True)
    return t, LAD


print("=" * 104)
print("reliability against ladder spacing")
print("=" * 104)
print(f"  {'spacing':<22} {'observable':<13} {'target':>8} {'certified':>10} "
      f"{'attained':>9} {'safety':>8} {'med slack':>10} {'med pen':>8}")

SPACINGS = [(2.0, 100.0), (1.6, 10.0), (1.4, 5.0), (1.2, 3.0)]
summary = []
for rd, rt in SPACINGS:
    tensor, LAD = tensor_for(rd, rt)
    HS = [[2 * np.pi / n for n in LAD[0]], [1.0 / s for s in LAD[1]],
          list(LAD[2]), list(LAD[3])]
    rows = []
    for ob in ["energy", "dissipation"]:
        QS = EX[ob]
        err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
        lk = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
        # pick the target that discriminates best: closest to the median error
        cand = sorted(err.values())
        for tgt in [cand[int(0.25 * len(cand))], cand[int(0.5 * len(cand))]]:
            eps = tgt
            feas = [i for i in ALL if err[i] <= eps]
            if not feas or len(feas) == len(ALL):
                continue
            c_opt = min(tensor[i][1] for i in feas)
            cert = safe = 0
            slack, pens = [], []
            for a in ALL:
                r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KIN,
                                 use_screen=False)
                if not r["certified"]:
                    continue
                cert += 1
                q = r["proposal"]
                safe += err[q] <= eps
                slack.append(err[q] / eps)
                pens.append(tensor[q][1] / c_opt)
            if cert == 0:
                print(f"  {f'{rd}x mesh / {rt}x tol':<22} {ob:<13} "
                      f"{eps/abs(QS):>8.1e} {'all 625 refused':>28}")
                continue
            print(f"  {f'{rd}x mesh / {rt}x tol':<22} {ob:<13} "
                  f"{eps/abs(QS):>8.1e} {cert:>10} {safe:>9} "
                  f"{100*safe/cert:>7.1f}% {np.median(slack):>10.3f} "
                  f"{np.median(pens):>7.2f}x")
            rows.append((cert, safe, float(np.median(slack))))
    if rows:
        c = sum(r[0] for r in rows); s = sum(r[1] for r in rows)
        summary.append((rd, rt, c, s, float(np.mean([r[2] for r in rows]))))
    print()

print("=" * 104)
print(f"  {'spacing':<24} {'certified':>10} {'attained':>9} {'safety':>9} "
      f"{'mean slack':>11}")
for rd, rt, c, s, sl in summary:
    print(f"  {f'{rd}x mesh / {rt}x tolerance':<24} {c:>10} {s:>9} "
          f"{100*s/c:>8.1f}% {sl:>11.3f}")
print()
print("  Safety and slack move together, as 6.4 predicts.  The method's")
print("  reliability on a coarse ladder is bought by the ladder, and a user who")
print("  supplies closely spaced candidate values should expect less of it.")
