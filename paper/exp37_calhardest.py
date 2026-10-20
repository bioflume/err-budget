"""If transfer fails only towards harder physics, calibrating at the hardest
member of the ensemble should make the whole ensemble safe.  That turns a
negative result into a usable rule, so it is worth testing rather than asserting.
"""
import os, pickle, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from errbudget import run_protocol
from errbudget.problems import Burgers1D

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
NAMES = ["N", "steps", "newton_tol", "gmres_tol"]
LAD = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
       [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
HS = [[2*np.pi/n for n in LAD[0]], [1.0/s for s in LAD[1]],
      list(LAD[2]), list(LAD[3])]
KIN = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
NUS = [0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.06,0.07,0.08,0.09,0.10,0.12,
       0.14,0.16,0.18,0.20,0.24,0.28,0.32,0.40,0.50]

# a tensor at the hardest member, so the method can be calibrated there
CACHE = os.path.join(D, "burgers_nu002_tensor.pkl")
import itertools, time
ALL = list(itertools.product(*[range(5) for _ in range(4)]))
p002 = Burgers1D(nu=0.02, T=1.0)
if os.path.exists(CACHE):
    t = pickle.load(open(CACHE, "rb"))
else:
    t, t0 = {}, time.time()
    for i, idx in enumerate(ALL):
        th = {NAMES[j]: LAD[j][v] for j, v in enumerate(idx)}
        obs, C, w, _ = p002.run(th)
        t[idx] = (obs, C)
        if i % 100 == 0: print(f"  {i}/625 {time.time()-t0:.0f}s", flush=True)
    pickle.dump(t, open(CACHE, "wb"))
    print(f"tensor at nu=0.02: {time.time()-t0:.0f}s", flush=True)

EX002 = p002.exact_observables()
print("="*92)
print("calibrate at the HARDEST ensemble member (nu = 0.02) and transport upward")
print("="*92)
for ob in ["energy", "dissipation"]:
    for tr in [1e-2, 3e-3]:
        QS = EX002[ob]
        look = lambda i, ob=ob: (t[tuple(i)][0][ob], t[tuple(i)][1])
        eps = tr*abs(QS)
        r = run_protocol(look, LAD, (2,2,3,3), eps, hs=HS, kinds=KIN, use_screen=False)
        rec = {NAMES[k]: LAD[k][j] for k, j in enumerate(r["proposal"])}
        nok, worst = 0, 0.0
        for nu in NUS:
            q = Burgers1D(nu=nu, T=1.0)
            ex = q.exact_observables()[ob]
            obs, c, _, _ = q.run(rec)
            rel = abs(obs[ob]-ex)/abs(ex)
            nok += rel <= tr; worst = max(worst, rel/tr)
        cost_here = t[r["proposal"]][1]
        print(f"  {ob:<12} target {tr:.0e}  ->  "
              + ", ".join(f"{k}={v:g}" for k,v in rec.items()))
        print(f"  {'':<12} transports to {nok}/{len(NUS)} members, worst margin "
              f"{worst:.2f}x, cost {cost_here:.3g} work units")
