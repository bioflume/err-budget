"""Experiments demanded by referee criticisms M4, M9c and M13.

M13  Is a reported residual of 5.8e-23 a measurement or round-off noise?
M9c  Does a recommendation calibrated at one physical parameter transfer to
     others?  The referee calls this the single most important missing
     experiment, because the payback argument depends on it.
M4   F_s = 1.25 was selected on the same problems on which safety is reported.
     Calibrate it on L0-L2 only and report L3/L4 as held out.
"""
import itertools
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from knobs import run_protocol
from knobs.problems import Burgers1D, Burgers2D, ReactionDiffusion2D

D = os.path.join(os.path.dirname(__file__), "data")

# ----------------------------------------------------------------- M13 ------
print("=" * 88)
print("M13  what is the resolution floor of our own observables?")
print("=" * 88)
b2 = Burgers2D(nu=0.05, T=1.0)
th = dict(N=64, steps=40, newton_tol=1e-8, gmres_tol=1e-8, ilu_drop=1e-3)
vals = [b2.run(th)[0]["gradient"] for _ in range(3)]
print(f"  same configuration, three evaluations: spread "
      f"{max(vals)-min(vals):.3e} (deterministic code, so this is 0)")
Q = vals[0]
eps_mach = np.finfo(float).eps
n = 64 * 64
floor_abs = abs(Q) * eps_mach * np.sqrt(n)
print(f"  observable magnitude {abs(Q):.4f}; float64 eps {eps_mach:.2e}; "
      f"{n} grid points")
print(f"  round-off floor of the summation ~ |Q| * eps * sqrt(N^2) = "
      f"{floor_abs:.3e} absolute, {floor_abs/abs(Q):.3e} relative")
print(f"  the residual reported for gmres_tol was 5.8e-23 absolute "
      f"= {5.8e-23/abs(Q):.1e} relative")
print(f"  -> that is {floor_abs/5.8e-23:.1e}x BELOW the round-off floor. It is a")
print("     fitted coefficient indistinguishable from zero, not a measurement.")
print("     Correct statement: the parameter's contribution is below the")
print(f"     resolution of the measurement, {floor_abs/abs(Q):.0e} relative.")

# ----------------------------------------------------------------- M9c ------
print()
print("=" * 88)
print("M9c  does a recommendation transfer across a PHYSICAL parameter?")
print("=" * 88)
LAD1 = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
        [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
NM1 = ["N", "steps", "newton_tol", "gmres_tol"]
HS1 = [[2 * np.pi / n for n in LAD1[0]], [1.0 / s for s in LAD1[1]],
       list(LAD1[2]), list(LAD1[3])]
KIN1 = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
t1 = pickle.load(open(os.path.join(D, "burgers_tensor.pkl"), "rb"))
ob = "dissipation"
b1 = Burgers1D(nu=0.05, T=1.0)
QS = b1.exact_observables()[ob]
look = lambda i: (t1[tuple(i)][0][ob], t1[tuple(i)][1])

target_rel = 3e-3
eps = target_rel * abs(QS)
start = (2, 2, 3, 3)
r = run_protocol(look, LAD1, start, eps, hs=HS1, kinds=KIN1, use_screen=False)
rec = {NM1[k]: LAD1[k][j] for k, j in enumerate(r["proposal"])}
print(f"  calibrated at nu = 0.05, relative target {target_rel:.0e}")
print(f"  recommendation: {rec}")
print()
print("  nu      exact Q      Q(recommended)   rel. error   target   verdict")
nus = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.15, 0.20]
n_ok = 0
for nu in nus:
    p = Burgers1D(nu=nu, T=1.0)
    ex = p.exact_observables()[ob]
    obs, cost, _, _ = p.run(rec)
    rel = abs(obs[ob] - ex) / abs(ex)
    ok = rel <= target_rel
    n_ok += ok
    print(f"  {nu:5.2f}  {ex:11.6f}  {obs[ob]:14.6f}  {rel:11.3e}  "
          f"{target_rel:.0e}   {'holds' if ok else 'FAILS'}")
print(f"\n  the recommendation transfers to {n_ok}/{len(nus)} physical settings")

# ----------------------------------------------------------------- M4 -------
print()
print("=" * 88)
print("M4  F_s calibrated on L2 only; L3 and L4 held out")
print("=" * 88)


def safety(tensor, LAD, HS, KIN, obs_list, exact, targets, fs, rel_target):
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    out = []
    for ob in obs_list:
        QS = exact[ob]
        lk = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
        tr = lambda i: abs(tensor[i][0][ob] - QS)
        for t in targets:
            e = t * abs(QS) if rel_target else t
            if not [i for i in ALL if tr(i) <= e]:
                continue
            c_opt = min(tensor[i][1] for i in ALL if tr(i) <= e)
            safe, pens = 0, []
            for a in ALL:
                rr = run_protocol(lk, LAD, a, e, hs=HS, kinds=KIN,
                                  use_screen=False, fs_cap=fs)
                q = rr["proposal"]
                safe += tr(q) <= e
                pens.append(tensor[q][1] / c_opt)
            out.append((ob, t, safe / len(ALL), float(np.median(pens)),
                        float(np.max(pens))))
    return out


LAD3 = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS3 = [[2 * np.pi / n for n in LAD3[0]], [1 / s for s in LAD3[1]],
       list(LAD3[2]), list(LAD3[3]), list(LAD3[4])]
KIN3 = ["discretisation", "discretisation", "nested-solver", "nested-solver",
        "cost-only"]
t3 = pickle.load(open(os.path.join(D, "burgers2d_tensor.pkl"), "rb"))

print("\n  F_s    L2 (calibration set)          L3 (held out)")
print("         safe   med pen  worst        safe   med pen  worst")
for fs in [1.0, 1.10, 1.25, 1.5, 2.0]:
    a = safety(t1, LAD1, HS1, KIN1, ["energy", "dissipation"],
               b1.exact_observables(), [1e-2, 3e-3, 1e-3], fs, False)
    b = safety(t3, LAD3, HS3, KIN3, ["energy", "gradient"],
               b2.exact_observables(), [1e-2, 3e-3], fs, True)
    sa = min(x[2] for x in a); pa = np.median([x[3] for x in a])
    wa = max(x[4] for x in a)
    sb = min(x[2] for x in b); pb = np.median([x[3] for x in b])
    wb = max(x[4] for x in b)
    print(f"  {fs:4.2f}   {100*sa:5.1f}%  {pa:6.2f}x  {wa:6.2f}x      "
          f"{100*sb:5.1f}%  {pb:6.2f}x  {wb:6.2f}x")
print("\n  (safe = worst cell; med pen = median over cells; worst = worst over all starts)")
