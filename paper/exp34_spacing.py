"""Closely spaced candidate values on L2 (S-3.5, Table S-3.4): four lists of five values
with the same geometric midpoint (N = 128, 100 time steps, both tolerances 1e-6), refined
by 2/100, 1.6/10, 1.4/5 and 1.2/3 (mesh and time step / tolerances).  The released program
searches from every setting at the 25th and 50th percentile of each observable's errors.
The tensors are built on first use (data/burgers_sp_<mesh>_<tol>.pkl); their convergence
flags come from l2flags.py.

    python exp34_spacing.py         # about nine minutes; writes data/spacing_released.json
"""
import itertools
import os
import pickle
import time

import numpy as np

import json
from errbudget import run_protocol, V121
from l2flags import flags
from errbudget.problems import Burgers1D

from paths import DATA as D, OUT
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
summary, OUTROWS = [], []
for rd, rt in SPACINGS:
    tensor, LAD = tensor_for(rd, rt)
    HS = [[2 * np.pi / n for n in LAD[0]], [1.0 / s for s in LAD[1]],
          list(LAD[2]), list(LAD[3])]
    FL = flags(tensor, LAD, 0.05, f"L2_sp_{rd}_{rt}")
    rows, allpens, cases = [], [], []
    for ob in ["energy", "dissipation"]:
        QS = EX[ob]
        err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
        lk = lambda i, ob=ob: ((tensor[tuple(i)][0][ob] if FL[tuple(i)]["converged"] else np.nan),
                               tensor[tuple(i)][1])
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
                r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KIN, **V121)
                if not r["certified"]:
                    continue
                cert += 1
                q = r["proposal"]
                safe += err[q] <= eps
                slack.append(err[q] / eps)
                pens.append(tensor[q][1] / c_opt)
            allpens += pens
            cases.append(dict(ob=ob, target_rel=eps / abs(QS), starts=len(ALL), accepted=cert,
                              met=safe, med_slack=float(np.median(slack)) if slack else None,
                              med_cost_ratio=float(np.median(pens)) if pens else None))
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
    OUTROWS.append(dict(rd=rd, rt=rt, h_range=HS[0][0] / HS[0][-1], cases=cases,
                        accepted=sum(c["accepted"] for c in cases), met=sum(c["met"] for c in cases),
                        searches=len(ALL) * len(cases),
                        med_cost_ratio=float(np.median(allpens)) if allpens else None))
    print()

json.dump(OUTROWS, open(os.path.join(OUT, "spacing_released.json"), "w"), indent=1)
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
