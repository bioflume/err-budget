"""Build the exhaustive tensors of the held-out problems H1 (scikit-fem) and H2
(OpenFOAM cavity, Re = 1000).  This collects solver output only; the program
is replayed on these tensors by exp62_heldout.py, with the configuration frozen
in PREREGISTRATION.md.  Resumable.
    python exp61_heldout_build.py H1|H2
"""
import itertools
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
which = sys.argv[1]
if which == "H1":
    import heldout_fem as M
    solver = M.run
else:
    import heldout_cavity as M
    _c = M.Cavity()
    solver = _c.run
CACHE = os.path.join(D, f"heldout_{which}.pkl")
part = CACHE + ".partial"
tensor = pickle.load(open(part, "rb")) if os.path.exists(part) else {}
ALL = list(itertools.product(*[range(len(l)) for l in M.LAD]))
t0 = time.time()
for n, idx in enumerate(ALL):
    if idx in tensor:
        continue
    th = {M.NAMES[k]: M.LAD[k][j] for k, j in enumerate(idx)}
    try:
        obs, cost, wall, info = solver(th)
    except Exception as exc:                       # a failed run is data too
        obs, cost, wall, info = ({k: float("nan") for k in M.EXACT if hasattr(M, "EXACT")}
                                 if hasattr(M, "EXACT") else {k: float("nan") for k in M.REF},
                                 float("inf"), 0.0, {"unavailable": str(exc)[:200]})
    tensor[idx] = (obs, cost, info)
    pickle.dump(tensor, open(part, "wb"))
    print(f"  {which} {len(tensor)}/{len(ALL)} {time.time()-t0:.0f}s {th} {info}", flush=True)
pickle.dump(tensor, open(CACHE, "wb"))
os.remove(part)
print(f"{which}: {len(tensor)} runs in {time.time()-t0:.0f}s")
