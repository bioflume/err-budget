"""Build the exhaustive OpenFOAM tensor for the DFG cylinder benchmark 2D-1.

7 mesh levels x 6 outer tolerances x 5 linear-solver tolerances x 4 GAMG
coarsest-level sizes = 840 simpleFoam runs.  The tensor is used only to grade
the method; the method itself sees only the runs it asks for.  Progress is
saved every few runs, so an interrupted build resumes where it stopped.
"""
import itertools
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
from openfoam_cylinder import DFGCylinder

NAMES = ["n_arc", "res_ctrl", "lin_tol", "n_coarsest"]
LAD = [[8, 12, 16, 24, 32, 48, 64],
       [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7],
       [1e-4, 1e-5, 1e-6, 1e-7, 1e-8],
       [10, 50, 250, 1000]]
ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
# v1.2: wall-cell centroids corrected and a second-order pressure difference
# added (openfoam_cylinder.py); the v1.1 tensor is kept as dfg_tensor.pkl
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                     os.environ.get("KNOBS_DFG_CACHE", "dfg_tensor_v12.pkl"))

if __name__ == "__main__":
    part = CACHE + ".partial"
    tensor = pickle.load(open(part, "rb")) if os.path.exists(part) else {}
    p = DFGCylinder()
    t0 = time.time()
    todo = [i for i in ALL if i not in tensor]
    print(f"{len(tensor)} runs cached, {len(todo)} to go", flush=True)
    for n, idx in enumerate(todo):
        th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
        obs, cost, wall, info = p.run(th)
        tensor[idx] = (obs, cost, info)
        if n % 10 == 0 or n == len(todo) - 1:
            pickle.dump(tensor, open(part, "wb"))
            print(f"  {len(tensor)}/{len(ALL)}  {time.time()-t0:.0f}s  last {th} "
                  f"{info.get('iters')} it {cost:.1f}s", flush=True)
    pickle.dump(tensor, open(CACHE, "wb"))
    os.remove(part)
    print(f"tensor: {len(tensor)} simpleFoam runs in {time.time()-t0:.0f}s")
