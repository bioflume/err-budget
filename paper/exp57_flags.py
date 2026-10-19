"""Convergence flags of the test-problem tensors (referee point 6 of revision 3).

Version 1.2 treats a run whose nested solver did not meet its own tolerance as
failed.  The stored tensors of L2-L5 predate the flags, so this script reruns
every configuration with the flagged solvers, checks that the observables and
the work-unit cost are bit-for-bit those of the stored tensor, and stores the
flags in data/flags_<problem>.pkl.  Resumable: partial results are kept.

    python exp57_flags.py L2|L3|L4|L5
"""
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grids
from knobs.problems import Burgers1D, Burgers2D, ReactionDiffusion2D, DrivenCavity

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
name = sys.argv[1]
P = getattr(grids, name)
prob = {"L2": lambda: Burgers1D(nu=0.05, T=1.0), "L3": lambda: Burgers2D(nu=0.05, T=1.0),
        "L4": ReactionDiffusion2D, "L5": DrivenCavity}[name]()
tensor = grids.load(P)
out_path = os.path.join(D, f"flags_{name}.pkl")
part = out_path + ".partial"
flags = pickle.load(open(part, "rb")) if os.path.exists(part) else {}
t0, mism = time.time(), 0
for n, idx in enumerate(grids.grid(P)):
    if idx in flags:
        continue
    th = {P["names"][k]: P["lad"][k][j] for k, j in enumerate(idx)}
    obs, cost, wall, info = prob.run(th)
    old = tensor[idx]
    same = all(obs[o] == old[0][o] for o in old[0]) and cost == old[1]
    mism += not same
    flags[idx] = dict(info, same=same)
    if n % 20 == 0:
        pickle.dump(flags, open(part, "wb"))
        print(f"  {name} {len(flags)}/{len(tensor)} {time.time()-t0:.0f}s, mismatches {mism}",
              flush=True)
pickle.dump(flags, open(out_path, "wb"))
if os.path.exists(part):
    os.remove(part)
bad = sum(not f["converged"] for f in flags.values())
print(f"{name}: {len(flags)} runs, {mism} differ from the stored tensor, "
      f"{bad} did not converge "
      f"(Newton: {sum(f.get('newton_fail', 0) > 0 or not f.get('newton_converged', True) for f in flags.values())}, "
      f"GMRES: {sum(f.get('gmres_fail', 0) > 0 for f in flags.values())})")
