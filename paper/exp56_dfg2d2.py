"""DFG benchmark 2D-2 (Re = 100) in OpenFOAM: analysis of the grid built by
exp56_dfg2d2_build.py (run `exp56_dfg2d2_build.py collect` first).

Observables: maximum drag and lift coefficients over a shedding period and the
Strouhal number, all of them extrema or a frequency rather than integrals.
The reference is the finest OpenFOAM run (n_arc = 64, dt = 2.5e-4); its
uncertainty is the larger change to the runs with twice the time step and with
the next coarser mesh, and a target is admitted only if it is at least ten
times that uncertainty.
"""
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from replay import grid, method_cell, baseline_cell, fmt_method
from exp47_transfer import failure_rate, PHIS
from errbudget import LocalBudget
import openfoam_dfg2d2 as F

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
z = pickle.load(open(os.path.join(D, "dfg2d2_tensor.pkl"), "rb"))
LAD, NAMES, ref = z["LAD"], z["NAMES"], z["ref"]
KINDS = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
HS = [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), list(LAD[3])]
OBS = ["cd_max", "cl_max", "st"]
ALLG = grid(LAD)
tensor = dict(z["tensor"])
missing = {i for i in ALLG if i not in tensor}
for i in missing:                       # not (yet) computed: treated as unusable
    tensor[i] = ({o: float("nan") for o in OBS}, float("inf"), {})
show = lambda i: "n_arc=%d, dt=%g, outer_tol=%.0e, p_tol=%.0e" % tuple(
    LAD[k][j] for k, j in enumerate(i))

print("=" * 100)
print(f"OpenFOAM pimpleFoam, DFG 2D-2 (Re = 100): {len(z['tensor'])}/{len(ALLG)} grid settings, "
      f"{len(ref)}/3 reference runs")
print("=" * 100)
if 0 not in ref:
    raise SystemExit("reference run not available yet")
R = ref[0][1]
unc = {}
for o in OBS:
    d = [abs(ref[k][1][o] - R[o]) / abs(R[o]) for k in (1, 2) if k in ref]
    unc[o] = max(d) if d else float("nan")
print("reference (n_arc 64, dt 2.5e-4):", {o: round(R[o], 6) for o in OBS})
print("relative uncertainty (dt x2, mesh one level coarser):", {o: f"{unc[o]:.1e}" for o in OBS})
print("published intervals (Schaefer & Turek 1996):", F.PUBLISHED)

usable = [i for i in ALLG if np.isfinite(tensor[i][1])]

# Is the reference in the asymptotic range?  Observed order and Richardson limit from
# n_arc 32 and 48 at dt = 5e-4 (tightest tolerances that ran) and the n_arc 64 run at
# dt = 5e-4 (reference run 1).
print("\nmesh convergence at dt = 5e-4 (n_arc 24, 32, 48, 64) and Richardson limit:")
from scipy.optimize import brentq
for o in OBS:
    seq = []
    for a in (2, 3, 4):
        i = None
        for c, d in [(3, 2), (3, 1), (2, 2), (2, 1)]:
            if np.isfinite(tensor[(a, 3, c, d)][1]):
                i = (a, 3, c, d); break
        seq.append((LAD[0][a], tensor[i][0][o] if i else np.nan))
    if 1 in ref:
        seq.append((64, ref[1][1][o]))
    line = f"  {o:<7} " + ", ".join(f"{n}: {q:.6f}" for n, q in seq)
    (n1, q1), (n2, q2), (n3, q3) = seq[-3:]
    r = (q2 - q1) / (q3 - q2)
    g = lambda p: (n1**-p - n2**-p) / (n2**-p - n3**-p) - r
    try:
        p = brentq(g, 0.05, 8.0)
        lim = q3 + (q3 - q2) * n3**-p / (n2**-p - n3**-p)
        line += f";  observed order {p:.2f}, limit {lim:.6f}, n_arc 64 off by {abs(q3-lim)/abs(lim):.1e}"
    except ValueError:
        line += ";  not monotone: no observed order"
    print(line)
print(f"usable settings: {len(usable)}; failed or not shedding: {len(ALLG)-len(usable)}")

# tolerance pairs from tightest to loosest; at large Courant numbers the tightest
# pair can fail, and the tightest pair that ran is used instead (marked *)
TOLS = sorted(((c, d) for c in range(len(LAD[2])) for d in range(len(LAD[3]))),
              key=lambda cd: (-(cd[0] + cd[1]), -cd[1]))


def tightest(a, b):
    for c, d in TOLS:
        if np.isfinite(tensor[(a, b, c, d)][1]):
            return (a, b, c, d)
    return None


print("\nrelative error at the tightest solver tolerances that ran (rows: n_arc; columns: dt;")
print("* = the tightest pair failed and a looser one is shown):")
for o in OBS:
    print(f"  {o}")
    for a in range(len(LAD[0])):
        row = []
        for b in range(len(LAD[1])):
            i = tightest(a, b)
            if i is None:
                miss = all((a, b, c, d) in missing for c, d in TOLS)
                row.append("    not run" if miss else "      fail ")
            else:
                mark = "*" if i[2:] != TOLS[0] else " "
                row.append(f"{abs(tensor[i][0][o]-R[o])/abs(R[o]):9.2e}{mark}")
        print(f"    n_arc {LAD[0][a]:>3}: " + " ".join(row))

print("\nfailed settings, by n_arc (rows) and dt (columns), out of 12 tolerance pairs:")
for a in range(len(LAD[0])):
    print(f"    n_arc {LAD[0][a]:>3}: " + " ".join(
        f"{sum(1 for c, d in TOLS if not np.isfinite(tensor[(a, b, c, d)][1]) and (a, b, c, d) not in missing):>3}"
        for b in range(len(LAD[1]))))

print("\npart of the error due to the solver tolerances (loosest against tightest that ran), by mesh, dt = 1e-3:")
for o in OBS:
    row = []
    for a in range(len(LAD[0])):
        i0, i1 = (a, 2, 0, 0), tightest(a, 2)
        if i1 is not None and np.isfinite(tensor[i0][1]):
            row.append(f"{abs(tensor[i0][0][o]-tensor[i1][0][o])/abs(R[o]):.1e}")
        else:
            row.append("fail")
    print(f"  {o:<7} " + "  ".join(row))

print("\nseparability of mesh and time step (Appendix: rho = |interaction| / larger main effect),")
print("adjacent levels, every usable tolerance pair:")
for o in OBS:
    rho = []
    for a in range(len(LAD[0]) - 1):
        for b in range(len(LAD[1]) - 1):
            for c, d in TOLS:
                q = [tensor[(a + x, b + y, c, d)] for x, y in ((0, 0), (1, 0), (0, 1), (1, 1))]
                if not all(np.isfinite(v[1]) for v in q):
                    continue
                v = [w[0][o] for w in q]
                I = (v[3] - v[2]) - (v[1] - v[0])
                M = max(abs(v[1] - v[0]), abs(v[2] - v[0]))
                if M > 1e-6 * abs(R[o]):
                    rho.append(abs(I) / M)
    rho = np.array(rho)
    print(f"  {o:<7} blocks {len(rho)}, median rho {np.median(rho):.2f}, "
          f"rho > 0.5 in {100*np.mean(rho > 0.5):.0f}%, rho > 1 in {100*np.mean(rho > 1):.0f}%")


def get(idx, ob):
    v = tensor[idx]
    return v[0][ob] if np.isfinite(v[1]) else np.nan
get.exact = R
print("\ntransfer of thresholds (as Table 3), pooled over the three observables (failure rate / comparisons):")
print(f"  {'parameter':<10} " + " ".join(f"phi={p:<6g}" for p in PHIS))
for k in range(4):
    row = []
    for phi in PHIS:
        f, n = failure_rate(get, LAD, k, OBS, phi, fallback=True)
        row.append(f"{100*f/max(n,1):>5.0f}% /{n:<3d}")
    print(f"  {NAMES[k]:<10} " + " ".join(row))

for o in OBS:
    for tg in [1e-2, 3e-3, 1e-3]:
        if not tg >= 2 * unc[o]:
            print(f"\n--- {o} {tg:.0e}: below twice the reference uncertainty; not graded")
            continue
        admitted = tg >= 10 * unc[o]
        eps = tg * abs(R[o])
        band = unc[o] / tg                   # reference uncertainty in units of eps
        print(f"\n--- {o}, target {tg:.0e} relative" + ("" if admitted else
              f"  [below ten times the reference uncertainty: graded with a band of "
              f"{band:.2f} eps]"))
        for lab, kw in [("method", {}), ("method, gated margin", dict(resid_margin=5.0))]:
            r = method_cell(tensor, LAD, HS, KINDS, o, R[o], eps, starts=usable, **kw)
            print(f"  {lab:<26} " + fmt_method(r))
            if not admitted:
                sure = sum(1 for m in r["misses"] if m[2] > 1 + band)
                print(f"      misses beyond the band (certain): {sure}; within it (possible): "
                      f"{len(r['misses']) - sure}")
            for m in r["misses"][:4]:
                print(f"      miss: {show(m[0])} -> {show(m[1])}: true {m[2]:.2f} eps, estimate {m[3]:.2f} eps")
        for lab, b in baseline_cell(tensor, LAD, o, R[o], eps, starts=usable).items():
            print(f"  {lab:<26} attained {b['attained']}/{b['n']}, penalty median "
                  f"{b['pen']:.2f}, runs median {b['runs']:.0f}")

print("\ngraded instead against FEATFLOW level 6 (", F.FEATFLOW, "):")
print("  the OpenFOAM reference differs from it by " + ", ".join(
    f"{o} {abs(R[o]-F.FEATFLOW[o])/F.FEATFLOW[o]:.1e}" for o in OBS))
for o in OBS:
    for tg in [3e-2, 1e-2]:
        Q = F.FEATFLOW[o]
        eps = tg * abs(Q)
        r = method_cell(tensor, LAD, HS, KINDS, o, Q, eps, starts=usable)
        b = baseline_cell(tensor, LAD, o, Q, eps, starts=usable)
        print(f"  {o:<7} {tg:.0e}: method " + fmt_method(r))
        print(f"  {'':<7}        sequential @eps attained {b['sequential @eps']['attained']}/"
              f"{b['sequential @eps']['n']}, coordinate search @eps/10 attained "
              f"{b['coordinate search @eps/10']['attained']}/{b['coordinate search @eps/10']['n']}")

# the misses of the drag at 3e-3: signs of the estimated mesh and time-step errors
from errbudget import run_protocol
eps = 3e-3 * abs(R["cd_max"])
lk = lambda i: (tensor[tuple(i)][0]["cd_max"], tensor[tuple(i)][1])
n_miss = n_opp = 0
for a in usable:
    r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KINDS)
    if r["certified"] and abs(tensor[r["proposal"]][0]["cd_max"] - R["cd_max"]) > eps:
        rk = r["budget"].rk
        n_miss += 1
        n_opp += bool(np.isfinite(rk[1]) and np.sign(rk[0]) != np.sign(rk[1]))
print(f"\ndrag at 3e-3: {n_miss} misses, {n_opp} of them with estimated mesh and time-step "
      "errors of opposite sign")

print("\nerror budget at a common user setting (n_arc 24, dt 2e-3, outer 1e-2, p 1e-6):")
for o in OBS:
    look = lambda i, o=o: (tensor[tuple(i)][0][o], tensor[tuple(i)][1])
    a = (2, 1, 0, 1)
    if not np.isfinite(tensor[a][1]):
        continue
    lb = LocalBudget(look, LAD, a, hs=HS)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    tot = sum(abs(x) for x in lb.rk if np.isfinite(x)) or 1
    print(f"  {o:<7} " + ", ".join(f"{NAMES[k]} {100*abs(lb.rk[k])/tot:.0f}%" for k in range(4))
          + f";  true error {abs(tensor[a][0][o]-R[o])/abs(R[o]):.2e}")
