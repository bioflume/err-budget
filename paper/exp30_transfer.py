"""Reuse of a setting across a parameter study on L2 (S-4.5, Table S-4.4): the program
searches at nu = 0.05 from N = 128, 100 time steps and both tolerances 1e-8 (released
program, unconverged runs failed), and the returned setting is run unchanged at 22
viscosities in [0.02, 0.5] and graded against the exact solution of each (Cole-Hopf,
evaluated in log space).

    python exp30_transfer.py        # about a minute; writes data/reuse_released.json
"""
import os
import pickle

import numpy as np

import json
from errbudget import run_protocol, V121
from errbudget.problems import Burgers1D

from paths import DATA as D, OUT
NAMES = ["N", "steps", "newton_tol", "gmres_tol"]
LAD = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
       [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
HS = [[2 * np.pi / n for n in LAD[0]], [1.0 / s for s in LAD[1]],
      list(LAD[2]), list(LAD[3])]
KIN = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
t1 = pickle.load(open(os.path.join(D, "burgers_tensor.pkl"), "rb"))
FL = pickle.load(open(os.path.join(D, "flags_L2.pkl"), "rb"))

NU_CAL = 0.05
NUS = [0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07, 0.08, 0.09,
       0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.24, 0.28, 0.32, 0.40, 0.50]

cal = Burgers1D(nu=NU_CAL, T=1.0)
EX = cal.exact_observables()

print("=" * 100)
print("transfer of a recommendation across a physical parameter (viscosity)")
print("=" * 100)
print(f"calibrated once at nu = {NU_CAL}; the SAME setting is then used at every")
print("other viscosity, with no further solver runs.")

RES = []
for ob, targets in [("energy", [1e-2, 3e-3]), ("dissipation", [1e-2, 3e-3])]:
    QS = EX[ob]
    look = lambda i, ob=ob: ((t1[tuple(i)][0][ob] if FL[tuple(i)]["converged"] else np.nan),
                             t1[tuple(i)][1])
    for tr in targets:
        eps = tr * abs(QS)
        start = (2, 2, 3, 3)
        r = run_protocol(look, LAD, start, eps, hs=HS, kinds=KIN, **V121)
        if not r["certified"]:
            print(f"\n--- {ob}, relative target {tr:.0e}: declined ({r.get('declined')})")
            RES.append(dict(ob=ob, target=tr, accepted=False, declined=r.get("declined")))
            continue
        rec = {NAMES[k]: LAD[k][j] for k, j in enumerate(r["proposal"])}
        print(f"\n--- {ob}, relative target {tr:.0e} "
              f"{'-'*(46-len(ob))}")
        print(f"  recommendation at nu={NU_CAL}: "
              + ", ".join(f"{k}={v:g}" for k, v in rec.items()))
        print(f"  {'nu':>6} {'exact':>12} {'transported':>13} {'rel. error':>11} "
              f"{'margin used':>12}   verdict")
        nok, worst, rows = 0, 0.0, []
        for nu in NUS:
            p = Burgers1D(nu=nu, T=1.0)
            ex = p.exact_observables()[ob]
            obs, cost, _, _ = p.run(rec)
            rel = abs(obs[ob] - ex) / abs(ex)
            ok = rel <= tr
            nok += ok
            worst = max(worst, rel / tr)
            rows.append(dict(nu=nu, rel=rel, ratio=rel / tr))
            print(f"  {nu:6.3f} {ex:12.6f} {obs[ob]:13.6f} {rel:11.3e} "
                  f"{rel/tr:11.2f}x   {'holds' if ok else 'FAILS'}")
        print(f"  -> the setting transports to {nok}/{len(NUS)} viscosities "
              f"spanning {min(NUS)}-{max(NUS)} ({max(NUS)/min(NUS):.0f}x)")
        RES.append(dict(ob=ob, target=tr, accepted=True, setting=rec, runs=r["runs"],
                        cost=t1[tuple(r["proposal"])][1], members=len(NUS), met=int(nok),
                        worst_ratio=worst, per_member=rows))

json.dump(RES, open(os.path.join(OUT, "reuse_released.json"), "w"), indent=1)
