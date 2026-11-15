"""Convergence flags for the stored L2 tensors that predate them (added for the reruns
with the released program, which treats a run whose nested solver missed its own
tolerance as failed; cf. exp57_flags.py, which did the same for the main L2-L5 grids).

flags(tensor, ladders, nu, tag) reruns every setting of the stored tensor, checks that
the observables and the work-unit cost are bit-for-bit those of the stored tensor, and
caches the flags in data/flags_<tag>.pkl.
"""
import itertools
import os
import pickle

from errbudget.problems import Burgers1D
from paths import OUT

NAMES = ["N", "steps", "newton_tol", "gmres_tol"]


def flags(tensor, LAD, nu, tag):
    path = os.path.join(OUT, f"flags_{tag}.pkl")
    if os.path.exists(path):
        return pickle.load(open(path, "rb"))
    p = Burgers1D(nu=nu, T=1.0)
    out, mism = {}, 0
    for idx in itertools.product(*[range(len(l)) for l in LAD]):
        th = {NAMES[j]: LAD[j][v] for j, v in enumerate(idx)}
        obs, cost, _, info = p.run(th)
        old = tensor[idx]
        same = all(obs[o] == old[0][o] for o in old[0]) and cost == old[1]
        mism += not same
        out[idx] = dict(info, same=same)
    pickle.dump(out, open(path, "wb"))
    print(f"flags {tag}: {len(out)} runs, {mism} differ from the stored tensor, "
          f"{sum(not f['converged'] for f in out.values())} did not converge", flush=True)
    return out
