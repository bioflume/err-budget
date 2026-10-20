"""The same baselines, on a production code rather than on our own problems."""
import itertools, os, pickle, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from errbudget import run_protocol
from baselines import seq, coord, ref_check
from lammps_madelung import LammpsMadelung

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
KINDS = ["discretisation", "discretisation", "cost-only"]

for ncell in [20, 40, 64]:
    path = os.path.join(D, f"lammps_tensor_n{ncell}.pkl")
    if not os.path.exists(path):
        print(f"n={ncell}: tensor not built"); continue
    tensor = pickle.load(open(path, "rb"))
    LAD = ([[6.0, 8.0, 10.0, 12.0, 14.0], [1e-2, 1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]]
           if ncell <= 20 else
           [[8.0, 10.0, 12.0, 14.0], [1e-3, 1e-4, 1e-5, 1e-6], [4, 5, 6, 7]])
    HS = [[1.0/c for c in LAD[0]], list(LAD[1]), [1.0/o for o in LAD[2]]]
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    p = LammpsMadelung(n_cells=ncell)
    QS = p.exact_observables()["energy_per_ion"]
    ok = [i for i in ALL if np.isfinite(tensor[i][0]["energy_per_ion"])]
    look = lambda i: (tensor[tuple(i)][0]["energy_per_ion"], tensor[tuple(i)][1])
    err = {i: abs(tensor[i][0]["energy_per_ion"] - QS)/abs(QS) for i in ok}

    print(f"\n{'='*96}")
    print(f"LAMMPS, {p.natoms:,} ions, {len(ok)}/{len(ALL)} configurations evaluable")
    print("=" * 96)
    for target in [1e-5, 1e-6]:
        feas = [i for i in ok if err[i] <= target]
        if not feas or len(feas) == len(ok): continue
        eps = target*abs(QS)
        c_opt = min(tensor[i][1] for i in feas)
        print(f"\n  target {target:.0e} relative")
        print(f"  {'procedure':<22} {'attained':>10} {'reliability':>12} "
              f"{'med penalty':>12} {'med runs':>9}")
        # method
        s=c=0; pens=[]; runs=[]
        for a in ok:
            r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, use_screen=False)
            if r["certified"]:
                c+=1; s += err[tuple(r["proposal"])] <= target
                pens.append(tensor[r["proposal"]][1]/c_opt)
            runs.append(r["runs"])
        print(f"  {'this method':<22} {s:>4}/{len(ok):<5} {100*s/len(ok):>11.1f}% "
              f"{np.median(pens) if pens else np.nan:>11.2f}x {np.median(runs):>9.0f}")
        for nm, fn in [("sequential", seq), ("coordinate search", coord)]:
            for frac, tag in [(1.0,"eps"), (0.1,"eps/10")]:
                s2=0; pn=[]; rn=[]
                for a in ok:
                    q, r = fn(look, LAD, a, eps*frac)
                    s2 += err.get(tuple(q), np.inf) <= target
                    pn.append(tensor[q][1]/c_opt); rn.append(r)
                print(f"  {nm+' @'+tag:<22} {s2:>4}/{len(ok):<5} "
                      f"{100*s2/len(ok):>11.1f}% {np.median(pn):>11.2f}x {np.median(rn):>9.0f}")
        q, r = ref_check(look, LAD, ok, eps)
        good = err.get(tuple(q), np.inf) <= target
        print(f"  {'reference acceptance':<22} {len(ok) if good else 0:>4}/{len(ok):<5} "
              f"{100.0*good:>11.1f}% {tensor[q][1]/c_opt:>11.2f}x {r:>9d}")
