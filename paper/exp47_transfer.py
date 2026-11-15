"""Failure rate of a convergence threshold chosen with every other parameter at its finest
value (helpers for exp64_transfer_floor.py).  The threshold of parameter k at a
combination c of the others is its coarsest value whose deviation from the value at k's
finest is within a fraction phi of the error remaining at c; a threshold chosen at the
finest combination fails at c when c needs a finer value of k.

    python exp47_transfer.py        # the rates without the adjustments of S-4.1
"""
import itertools
import os
import pickle

import numpy as np

import grids

PHIS = [0.01, 0.05, 0.1, 0.25, 0.5]


def failure_rate(get, LAD, k, obs, phi, fallback=False):
    """With fallback=True, an observable whose sweep at the finest setting of the other
    parameters contains a failed run is calibrated instead at the most refined setting
    of the other parameters whose sweep has none (used where fine meshes fail at large
    time steps)."""
    others = [j for j in range(len(LAD)) if j != k]
    combos = list(itertools.product(*[range(len(LAD[j])) for j in others]))
    fin = tuple(len(LAD[j]) - 1 for j in others)

    def thr(c, ob):
        def at(j):
            idx = [0] * len(LAD)
            for kk, v in zip(others, c):
                idx[kk] = v
            idx[k] = j
            return get(tuple(idx), ob)
        vals = [at(j) for j in range(len(LAD[k]))]
        if not all(np.isfinite(v) for v in vals):
            return None
        fine = vals[-1]
        tot = abs(fine - get.exact[ob])
        if tot <= 0:
            return 0
        for j, v in enumerate(vals):
            if abs(v - fine) <= phi * tot:
                return j
        return len(LAD[k]) - 1

    n = fail = 0
    for ob in obs:
        ref = thr(fin, ob)
        if ref is None and fallback:
            for c in sorted(combos, key=lambda c: (-sum(c), [-x for x in c])):
                ref = thr(c, ob)
                if ref is not None:
                    break
        if ref is None:
            continue
        for c in combos:
            t = thr(c, ob)
            if t is None:
                continue
            n += 1
            fail += t > ref
    return fail, n


def table(name, get, LAD, NAMES, KINDS, obs):
    print(f"\n{name}: failure rate of a threshold calibrated at the finest setting "
          f"of the other parameters")
    print(f"  {'parameter':<14} {'class':<15} " + " ".join(f"phi={p:<6g}" for p in PHIS))
    for k in range(len(LAD)):
        row = []
        for phi in PHIS:
            f, n = failure_rate(get, LAD, k, obs, phi)
            row.append(f"{100*f/max(n,1):>9.0f}%")
        print(f"  {NAMES[k]:<14} {KINDS[k]:<15} " + " ".join(row))


def from_tensor(t, exact):
    def get(idx, ob):
        v = t[idx]
        return v[0][ob] if np.isfinite(v[1]) else np.nan
    get.exact = exact
    return get


if __name__ == "__main__":
    D = os.path.join(os.path.dirname(__file__), "data")
    probs = [grids.L2, grids.L3, grids.L4]
    if os.path.exists(os.path.join(D, grids.L5["file"])):
        probs.append(grids.L5)
    for P in probs:
        t = grids.load(P)
        table(P["name"], from_tensor(t, grids.exact(P)), P["lad"], P["names"], P["kinds"], P["obs"])

    if os.path.exists(os.path.join(D, "dfg_tensor.pkl")):
        from exp46_dfg_build import LAD as DLAD, NAMES as DN
        from openfoam_cylinder import REF
        t = pickle.load(open(os.path.join(D, "dfg_tensor.pkl"), "rb"))
        table("DFG 2D-1 (OpenFOAM)", from_tensor(t, REF), DLAD, DN,
              ["discretisation", "nested-solver", "nested-solver", "cost-only"],
              ["cd", "dp"])
