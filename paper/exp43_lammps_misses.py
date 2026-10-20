"""Section 8.1: the certified misses on LAMMPS at 262,144 ions, and why."""
import sys, os, pickle, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from errbudget import run_protocol
from lammps_madelung import LammpsMadelung

N = 64
t = pickle.load(open(os.path.join(os.path.dirname(__file__), "data", f"lammps_tensor_n{N}.pkl"), "rb"))
QS = LammpsMadelung(n_cells=N).exact_observables()["energy_per_ion"]
LAD = [[8.0, 10.0, 12.0, 14.0], [1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]]
HS = [[1/c for c in LAD[0]], list(LAD[1]), [1/o for o in LAD[2]]]
KINDS = ["discretisation", "discretisation", "cost-only"]
look = lambda i: (t[tuple(i)][0]["energy_per_ion"], t[tuple(i)][1])
ok = [i for i in t if np.isfinite(t[i][1])]
show = lambda i: "%g A, %.0e, order %d" % tuple(LAD[k][j] for k, j in enumerate(i))
for tg in [1e-5, 1e-6]:
    eps = tg * abs(QS); cert = miss = 0
    for a in ok:
        r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, use_screen=False)
        if not r["certified"]:
            continue
        cert += 1; q = r["proposal"]
        e = abs(t[q][0]["energy_per_ion"] - QS)
        if e > eps:
            miss += 1
            fin = []
            for k in range(3):
                j = q[k] + 1
                if j >= len(LAD[k]):
                    fin.append(f"{['cutoff','accuracy','order'][k]}: end of ladder")
                else:
                    idx = list(q); idx[k] = j
                    fin.append(f"{['cutoff','accuracy','order'][k]}: finer setting "
                               + ("unusable" if not np.isfinite(t[tuple(idx)][1]) else "usable"))
            print(f"  target {tg:.0e}: start {show(a)} -> certified {show(q)}, true error "
                  f"{e/eps:.2f} eps, estimate {abs(r['budget'].E0)/eps:.2f} eps; " + "; ".join(fin))
    print(f"target {tg:.0e}: {len(ok)} starts, {cert} certified, {miss} missed")
