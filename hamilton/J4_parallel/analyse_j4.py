"""Analysis of J4: cost ordering and recommendations against the number of ranks.
    python analyse_j4.py [config]
"""
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
import parallel_hpc as J
import exp55_parallel as E
from lammps_madelung import MADELUNG, KE
from knobs import run_protocol
import pooled as P

RES = os.path.join(os.environ.get("WORK", "."), "parallel_hpc")
CFG = P.CONFIGS[sys.argv[1] if len(sys.argv) > 1 else "v1.2"]


def load(prefix):
    out = {}
    for R in J.RANKS:
        f = os.path.join(RES, f"{prefix}_R{R}", "tensor.pkl")
        if os.path.exists(f):
            out[R] = pickle.load(open(f, "rb"))
    return out


def compare(name, tensors, LAD, HS, KINDS, ob, QS, targets):
    keys = sorted(tensors)
    if len(keys) < 2:
        print(f"{name}: fewer than two widths finished"); return
    ALL = [i for i in tensors[keys[0]] if all(np.isfinite(tensors[k][i][1]) for k in keys)]
    print(f"\n{name}: {len(ALL)} settings usable at every width {keys}")
    base = tensors[keys[0]]
    rk = lambda v: np.argsort(np.argsort(v))
    for k in keys:
        t = tensors[k]
        c0 = np.log([base[i][1] for i in ALL]); c1 = np.log([t[i][1] for i in ALL])
        rho = np.corrcoef(rk(c0), rk(c1))[0, 1]
        print(f"  {k:>3} ranks: median cost {np.median([t[i][1] for i in ALL])*1e3:9.2f} ms, "
              f"rank correlation of costs with {keys[0]} rank(s): {rho:.3f}")
    for tg in targets:
        eps = tg * abs(QS)
        recs = {}
        for k in keys:
            t = tensors[k]
            lk = lambda i, t=t: (t[tuple(i)][0][ob], t[tuple(i)][1])
            recs[k] = {a: run_protocol(lk, LAD, a, eps, hs=HS, kinds=KINDS, **CFG["kw"]) for a in ALL}
            acc = [r for r in recs[k].values() if r["certified"]]
            miss = sum(abs(t[tuple(r["proposal"])][0][ob] - QS) > eps for r in acc)
            feas = [i for i in ALL if abs(t[i][0][ob] - QS) <= eps]
            opt = min((t[i][1] for i in feas), default=np.nan)
            pen = [t[tuple(r["proposal"])][1] / opt for r in acc]
            print(f"  target {tg:.0e}, {k:>3} ranks: accepted {len(acc)}/{len(ALL)}, missed {miss}, "
                  f"median penalty {np.median(pen) if pen else np.nan:.2f}")
        for k in keys[1:]:
            same = sum(tuple(recs[k][a]["proposal"]) == tuple(recs[keys[0]][a]["proposal"])
                       and recs[k][a]["certified"] == recs[keys[0]][a]["certified"] for a in ALL)
            print(f"    same outcome at {k} ranks as at {keys[0]}: {same}/{len(ALL)}")


QL = -0.5 * MADELUNG * KE / E.R0
compare("LAMMPS, 262,144 ions", load("lmp"), E.LADL, E.HSL, E.KL, "energy_per_ion", QL, [1e-4, 1e-5])
fr = os.path.join(RES, "gmx_ref", "ref.json")
if os.path.exists(fr):
    ref = json.load(open(fr))
    print(f"\nGROMACS reference {ref['QS']:.6f} kJ/mol (spread of two converged settings "
          f"{ref['spread']:.2e})")
    HSG = [[1.0 / c for c in J.LADG[0]], list(J.LADG[1]), [1.0 / o for o in J.LADG[2]], list(J.LADG[3])]
    compare(f"GROMACS double, {216*J.NBOX**3} waters", load("gmx"), J.LADG, HSG,
            ["discretisation", "discretisation", "cost-only", "discretisation"], "potential",
            ref["QS"], [1e-5, 1e-6])
