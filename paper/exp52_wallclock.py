"""Item: does the work-unit cost model matter?  Table 5 with measured time.

L2 and L3 count cost in work units (GMRES iterations times nonzeros).  Here
every configuration is re-run and timed by the CPU time of the process
(time.process_time, insensitive to other load on the machine), and the method
is replayed from every start with that time as the cost.  Reported per cell:
how many starts end at the same recommendation under both cost models, and the
median penalty when each is graded by its own cost model and by CPU time.
"""
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import grids
from errbudget import run_protocol
from errbudget.problems import Burgers1D, Burgers2D

D = os.path.join(os.path.dirname(__file__), "data")
MAKERS = {"L2": lambda: Burgers1D(nu=0.05, T=1.0), "L3": lambda: Burgers2D(nu=0.05, T=1.0)}

for P in [grids.L2, grids.L3]:
    t = grids.load(P)
    path = os.path.join(D, f"cputime_{P['name']}.pkl")
    if os.path.exists(path):
        cpu = pickle.load(open(path, "rb"))
    else:
        prob, cpu, t0 = MAKERS[P["name"]](), {}, time.time()
        for idx in grids.grid(P):
            th = {P["names"][k]: P["lad"][k][j] for k, j in enumerate(idx)}
            c0 = time.process_time()
            obs, C, w, _ = prob.run(th)
            cpu[idx] = time.process_time() - c0
            assert abs(obs[P["obs"][0]] - t[idx][0][P["obs"][0]]) <= 1e-9 * abs(t[idx][0][P["obs"][0]]) + 1e-14
        pickle.dump(cpu, open(path, "wb"))
        print(f"{P['name']}: {len(cpu)} runs re-timed in {time.time()-t0:.0f}s")
    EX = grids.exact(P)
    HS = grids.hs(P)
    ALL = grids.grid(P)
    rho = np.corrcoef(np.log([t[i][1] for i in ALL]), np.log([cpu[i] for i in ALL]))[0, 1]
    print(f"\n{P['name']}: correlation of log work units with log CPU time over the grid: {rho:.3f}")
    print(f"  {'cell':<26} {'same rec.':>10} {'cert (wu/cpu)':>14} {'missed':>7} "
          f"{'pen wu->wu':>11} {'pen cpu->cpu':>12} {'pen wu->cpu':>12}")
    for ob in P["obs"]:
        for tg in P["targets"]:
            eps = tg * abs(EX[ob]) if P["rel"] else tg
            err = {i: abs(t[i][0][ob] - EX[ob]) for i in ALL}
            feas = [i for i in ALL if err[i] <= eps]
            if not feas:
                continue
            opt_wu = min(t[i][1] for i in feas); opt_cpu = min(cpu[i] for i in feas)
            lk_wu = lambda i: (t[tuple(i)][0][ob], t[tuple(i)][1])
            lk_cpu = lambda i: (t[tuple(i)][0][ob], cpu[tuple(i)])
            same = cw = cc = miss = 0; p1, p2, p3 = [], [], []
            for a in ALL:
                r1 = run_protocol(lk_wu, P["lad"], a, eps, hs=HS, kinds=P["kinds"])
                r2 = run_protocol(lk_cpu, P["lad"], a, eps, hs=HS, kinds=P["kinds"])
                same += r1["proposal"] == r2["proposal"] and r1["certified"] == r2["certified"]
                if r1["certified"]:
                    cw += 1; p1.append(t[r1["proposal"]][1] / opt_wu)
                    p3.append(cpu[r1["proposal"]] / opt_cpu)
                if r2["certified"]:
                    cc += 1; p2.append(cpu[r2["proposal"]] / opt_cpu)
                    miss += err[r2["proposal"]] > eps
            print(f"  {ob + ' ' + format(tg, '.0e'):<26} {same:>5}/{len(ALL):<4} "
                  f"{cw:>6}/{cc:<6} {miss:>7} {np.median(p1):>11.2f} {np.median(p2):>12.2f} "
                  f"{np.median(p3):>12.2f}")
