"""How far the additive error model holds, and whether the estimate covers its failure
(S-3.4).  For every 2x2 design of adjacent values of two parameters k != l, at a fixed
setting z of the others,

    I = [E(b,d;z) - E(a,d;z)] - [E(b,c;z) - E(a,c;z)],
    M = max(|E(b,c;z) - E(a,c;z)|, |E(a,d;z) - E(a,c;z)|),   rho = |I| / M,

over the designs with M above a tenth of the tightest graded target, on L2-L4.  Also the
ratio of the estimate to the true error at every setting, against the local rho.

    python exp28_separability.py
"""
import itertools
import os
import pickle

import numpy as np

from errbudget import LocalBudget
from errbudget.problems import Burgers1D, Burgers2D, ReactionDiffusion2D

D = os.path.join(os.path.dirname(__file__), "data")
FLOOR = 1e-14          # nothing below the double-precision resolution counts


def rho_distribution(tensor, LAD, ob, QS, seed=0):
    """All 2x2 second differences of the error tensor, normalised."""
    n = len(LAD)
    E = {i: tensor[i][0][ob] - QS for i in tensor}
    scale = abs(QS) * FLOOR
    rng = np.random.default_rng(seed)
    out = []
    pairs = list(itertools.combinations(range(n), 2))
    rest = lambda k, l: [m for m in range(n) if m not in (k, l)]
    for k, l in pairs:
        for a, b in itertools.combinations(range(len(LAD[k])), 2):
            for c, d in itertools.combinations(range(len(LAD[l])), 2):
                zs = list(itertools.product(*[range(len(LAD[m]))
                                              for m in rest(k, l)]))
                if len(zs) > 40:
                    zs = [zs[i] for i in rng.choice(len(zs), 40, replace=False)]
                for z in zs:
                    def idx(ik, il):
                        v = [0] * n
                        v[k], v[l] = ik, il
                        for m, val in zip(rest(k, l), z):
                            v[m] = val
                        return tuple(v)
                    try:
                        Eac, Ebc = E[idx(a, c)], E[idx(b, c)]
                        Ead, Ebd = E[idx(a, d)], E[idx(b, d)]
                    except KeyError:
                        continue
                    if not all(np.isfinite(x) for x in (Eac, Ebc, Ead, Ebd)):
                        continue
                    I = (Ebd - Ead) - (Ebc - Eac)
                    M = max(abs(Ebc - Eac), abs(Ead - Eac))
                    if M <= scale:
                        continue
                    out.append(abs(I) / M)
    return np.array(out)


def estimate_quality(tensor, LAD, HS, ob, QS, fs_cap=1.25):
    """Is the estimate conservative, and does it fail where rho is large?"""
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    look = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    rows = []
    for a in ALL:
        Etrue = abs(tensor[a][0][ob] - QS)
        if Etrue <= abs(QS) * FLOOR:
            continue
        lb = LocalBudget(look, LAD, a, hs=HS, fs_cap=fs_cap)
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        if not np.isfinite(lb.E0):
            continue
        rows.append((abs(lb.E0) / Etrue, abs(lb.E0) / max(lb.Fs, 1e-9) / Etrue,
                     local_rho(tensor, LAD, ob, QS, a)))
    return np.array(rows)


def local_rho(tensor, LAD, ob, QS, a):
    """The worst 2x2 interaction in the anchor's immediate neighbourhood."""
    n = len(LAD)
    E = lambda i: tensor[i][0][ob] - QS
    scale = abs(QS) * FLOOR
    worst = 0.0
    for k, l in itertools.combinations(range(n), 2):
        if a[k] + 1 >= len(LAD[k]) or a[l] + 1 >= len(LAD[l]):
            continue
        def idx(dk, dl):
            v = list(a); v[k] += dk; v[l] += dl
            return tuple(v)
        try:
            Eac, Ebc, Ead, Ebd = (E(idx(0, 0)), E(idx(1, 0)),
                                  E(idx(0, 1)), E(idx(1, 1)))
        except KeyError:
            continue
        if not all(np.isfinite(x) for x in (Eac, Ebc, Ead, Ebd)):
            continue
        I = (Ebd - Ead) - (Ebc - Eac)
        M = max(abs(Ebc - Eac), abs(Ead - Eac))
        if M > scale:
            worst = max(worst, abs(I) / M)
    return worst


def q(x, p):
    return float(np.percentile(x, p)) if len(x) else float("nan")


print("=" * 100)
print("(1) the separability residual rho, defined and distributed")
print("=" * 100)
print(f"{'problem / observable':<34} {'n':>7} {'median':>8} {'p90':>8} "
      f"{'p99':>8} {'max':>9}  {'frac>1':>7}")

b1 = Burgers1D(nu=0.05, T=1.0)
t1 = pickle.load(open(os.path.join(D, "burgers_tensor.pkl"), "rb"))
LAD1 = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
        [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
HS1 = [[2 * np.pi / n for n in LAD1[0]], [1.0 / s for s in LAD1[1]],
       list(LAD1[2]), list(LAD1[3])]
E1 = b1.exact_observables()

b2 = Burgers2D(nu=0.05, T=1.0)
t2 = pickle.load(open(os.path.join(D, "burgers2d_tensor.pkl"), "rb"))
LAD2 = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS2 = [[2 * np.pi / n for n in LAD2[0]], [1.0 / s for s in LAD2[1]],
       list(LAD2[2]), list(LAD2[3]), list(LAD2[4])]
E2 = b2.exact_observables()

p4 = ReactionDiffusion2D()
t4 = pickle.load(open(os.path.join(D, "reactdiff_tensor.pkl"), "rb"))
LAD4 = [[16, 32, 64, 128], [5, 10, 20, 40], [1, 2, 4, 8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS4 = [[2 * np.pi / n for n in LAD4[0]], [1.0 / s for s in LAD4[1]],
       [1.0 / r for r in LAD4[2]], list(LAD4[3]), list(LAD4[4])]
E4 = p4.exact_observables()

CASES = [("L2 1D Burgers / energy", t1, LAD1, HS1, "energy", E1),
         ("L2 1D Burgers / dissipation", t1, LAD1, HS1, "dissipation", E1),
         ("L3 2D Burgers / energy", t2, LAD2, HS2, "energy", E2),
         ("L3 2D Burgers / gradient", t2, LAD2, HS2, "gradient", E2),
         ("L4 Gray-Scott / inhibitor_mass", t4, LAD4, HS4, "inhibitor_mass", E4),
         ("L4 Gray-Scott / pattern_sharpness", t4, LAD4, HS4,
          "pattern_sharpness", E4)]

RHO = {}
for name, T, LAD, HS, ob, EX in CASES:
    r = rho_distribution(T, LAD, ob, EX[ob])
    RHO[name] = r
    print(f"{name:<34} {len(r):>7} {q(r,50):>8.3f} {q(r,90):>8.3f} "
          f"{q(r,99):>8.3f} {r.max() if len(r) else np.nan:>9.2f}  "
          f"{100*np.mean(r>1):>6.1f}%")

print("\n  rho = 0 would be exact additivity.  The model is a search heuristic,")
print("  and these numbers say how badly it is violated: the cross term reaches")
print("  the size of the main effect and beyond in a minority of blocks.")

print()
print("=" * 100)
print("(2) is the estimate conservative, and does it fail where rho is large?")
print("=" * 100)
print(f"{'problem / observable':<34} {'anchors':>8} {'Ehat/Etrue':>22} "
      f"{'unmargined <1':>14} {'margined <1':>12}")
print(f"{'':<34} {'':>8} {'p10 / median / p90':>22}")
for name, T, LAD, HS, ob, EX in CASES:
    r = estimate_quality(T, LAD, HS, ob, EX[ob])
    if not len(r):
        continue
    ratio, raw, rho = r[:, 0], r[:, 1], r[:, 2]
    print(f"{name:<34} {len(r):>8} "
          f"{q(ratio,10):>6.2f} /{q(ratio,50):>6.2f} /{q(ratio,90):>7.2f} "
          f"{100*np.mean(raw<1):>13.1f}% {100*np.mean(ratio<1):>11.1f}%")
    # cross-tabulate against the local interaction strength
    lo, hi = rho <= 0.1, rho > 0.5
    if lo.sum() and hi.sum():
        print(f"{'    rho <= 0.1 / rho > 0.5':<34} "
              f"{lo.sum():>4}/{hi.sum():<4} "
              f"median Ehat/Etrue {np.median(ratio[lo]):>6.2f} / "
              f"{np.median(ratio[hi]):.2f};  under-estimates "
              f"{100*np.mean(ratio[lo]<1):.1f}% / {100*np.mean(ratio[hi]<1):.1f}%")

print()
print("=" * 100)
print("(3) the non-transferability mechanism: does e_tol/tol grow with dt?")
print("=" * 100)
print("  1D Burgers at N = 128, dissipation.  e_tol(dt, tol) is the deviation")
print("  caused by loosening gmres_tol from 1e-10 to tol, at fixed dt.")
print()
kN = LAD1[0].index(128)
kn = LAD1[2].index(1e-8)
fine = len(LAD1[3]) - 1
print(f"  {'steps':>6} {'dt':>8} {'matvecs':>9}  " +
      "  ".join(f"e/tol @{t:.0e}" for t in LAD1[3][:3]))
for js, s in enumerate(LAD1[1]):
    base = t1[(kN, js, kn, fine)][0]["dissipation"]
    cells = []
    for jt in range(3):
        tol = LAD1[3][jt]
        v = t1[(kN, js, kn, jt)][0]["dissipation"]
        cells.append(abs(v - base) / tol)
    _, _, _, info = b1.run(dict(N=128, steps=s, newton_tol=1e-8, gmres_tol=1e-8))
    print(f"  {s:>6} {1.0/s:>8.4f} {info['matvecs']:>9}  " +
          "  ".join(f"{c:>11.3e}" for c in cells))
print()
print("  If the ratios rise down each column, the tolerance's contribution to")
print("  the observable grows as the step coarsens -- which is what makes a")
print("  threshold calibrated at one discretisation invalid at another.")
