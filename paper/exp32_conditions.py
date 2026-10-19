"""Under WHICH conditions does a calibrated tolerance threshold fail to transfer?

The earlier claim was a law -- that nested-solver influence does not transfer
across discretisation settings.  The evidence supports a frequency, and a
frequency is only useful with its conditions attached.  So: compute, for every
combination of the other parameters, the masking threshold of the GMRES
tolerance, and cross-tabulate the failures of a threshold calibrated at the
finest combination against what distinguishes those combinations.

The threshold of parameter k at a configuration is the coarsest value of k whose
deviation from the fully refined value of k stays within a fraction of the total
error there.  A threshold "transfers" to another configuration if it is still at
least as fine as that configuration's own threshold; it fails if it is too coarse,
which is the unsafe direction.
"""
import itertools
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from knobs.problems import Burgers1D

D = os.path.join(os.path.dirname(__file__), "data")
LAD = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
       [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
NAMES = ["N", "steps", "newton_tol", "gmres_tol"]
t1 = pickle.load(open(os.path.join(D, "burgers_tensor.pkl"), "rb"))
b1 = Burgers1D(nu=0.05, T=1.0)
EX = b1.exact_observables()

K = 3                      # the parameter whose threshold we calibrate: gmres_tol
OTHERS = [0, 1, 2]
FINEST = tuple(len(LAD[k]) - 1 for k in OTHERS)


def threshold(combo, ob, frac=0.1):
    """Coarsest gmres_tol whose deviation stays under frac of the total error."""
    QS = EX[ob]
    def at(j):
        idx = [0, 0, 0, 0]
        for k, v in zip(OTHERS, combo):
            idx[k] = v
        idx[K] = j
        return t1[tuple(idx)][0][ob]
    fine = at(len(LAD[K]) - 1)
    tot = abs(fine - QS)
    if tot <= 0:
        return 0
    for j in range(len(LAD[K])):
        if abs(at(j) - fine) <= frac * tot:
            return j
    return len(LAD[K]) - 1


print("=" * 98)
print("a threshold calibrated with everything else at its finest, applied "
      "everywhere else")
print("=" * 98)
combos = list(itertools.product(*[range(len(LAD[k])) for k in OTHERS]))
for ob in ["energy", "dissipation"]:
    ref = threshold(FINEST, ob)
    # ladder index rises with tightness, so the calibrated value is TOO
    # COARSE exactly when some other combination needs a higher index.
    fails = [c for c in combos if threshold(c, ob) > ref]
    print(f"\n  {ob}: threshold at the finest combination is "
          f"gmres_tol = {LAD[K][ref]:.0e}")
    print(f"    too coarse in {len(fails)}/{len(combos)} of the other "
          f"combinations ({100*len(fails)/len(combos):.0f}%)")
    print(f"    failure rate broken down by the condition that distinguishes "
          f"them:")
    for ax, nm in [(0, "N"), (1, "steps"), (2, "newton_tol")]:
        print(f"      by {nm}:")
        for j, v in enumerate(LAD[ax]):
            sub = [c for c in combos if c[OTHERS.index(ax)] == j]
            f = [c for c in sub if threshold(c, ob) > ref]
            print(f"        {nm} = {v:<8g}  {len(f):>3}/{len(sub):<3} fail "
                  f"({100*len(f)/len(sub):>3.0f}%)")

print()
print("=" * 98)
print("the same question for the classes that are supposed to be safe")
print("=" * 98)


def threshold_k(combo, k, ob, frac=0.1):
    QS = EX[ob]
    others = [j for j in range(4) if j != k]
    def at(j):
        idx = [0, 0, 0, 0]
        for kk, v in zip(others, combo):
            idx[kk] = v
        idx[k] = j
        return t1[tuple(idx)][0][ob]
    fine = at(len(LAD[k]) - 1)
    tot = abs(fine - QS)
    if tot <= 0:
        return 0
    for j in range(len(LAD[k])):
        if abs(at(j) - fine) <= frac * tot:
            return j
    return len(LAD[k]) - 1


print(f"  {'parameter':<14} {'class':<16} {'observable':<12} "
      f"{'too coarse elsewhere':>22}")
for k, nm, cls in [(0, "N", "discretisation"), (1, "steps", "discretisation"),
                   (2, "newton_tol", "nested solver"),
                   (3, "gmres_tol", "nested solver")]:
    others = [j for j in range(4) if j != k]
    cmb = list(itertools.product(*[range(len(LAD[j])) for j in others]))
    fin = tuple(len(LAD[j]) - 1 for j in others)
    for ob in ["energy", "dissipation"]:
        ref = threshold_k(fin, k, ob)
        f = [c for c in cmb if threshold_k(c, k, ob) > ref]
        print(f"  {nm:<14} {cls:<16} {ob:<12} "
              f"{len(f):>10}/{len(cmb):<4} ({100*len(f)/len(cmb):>3.0f}%)")

print()
print("=" * 98)
print("Table 3: failure rate pooled over both observables, across the fraction phi")
print("=" * 98)
PHIS = [0.01, 0.05, 0.1, 0.25, 0.5]
print(f"  {'parameter':<14} " + " ".join(f"phi={p:<5g}" for p in PHIS))
for k, nm in [(0, "N"), (1, "steps"), (2, "newton_tol"), (3, "gmres_tol")]:
    others = [j for j in range(4) if j != k]
    cmb = list(itertools.product(*[range(len(LAD[j])) for j in others]))
    fin = tuple(len(LAD[j]) - 1 for j in others)
    row = []
    for frac in PHIS:
        n = 0
        for ob in ["energy", "dissipation"]:
            ref = threshold_k(fin, k, ob, frac)
            n += sum(threshold_k(c, k, ob, frac) > ref for c in cmb)
        row.append(f"{100*n/(2*len(cmb)):>8.0f}%")
    print(f"  {nm:<14} " + " ".join(row))
