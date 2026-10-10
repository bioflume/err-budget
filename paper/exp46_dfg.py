"""Item 1: the method on OpenFOAM, DFG benchmark 2D-1 (drag, lift, pressure drop).

Reads the tensor built by exp46_dfg_build.py and reports
  * the error of each observable against mesh level at the tightest tolerances;
  * the error budget at the OpenFOAM tutorial settings, and the true split;
  * the transfer of thresholds (Table 3) on this code;
  * the method from every start, with the baselines;
  * the recommendation from the tutorial settings against the settings themselves;
  * a falsified declaration: relax_U declared cost-only (live runs).
"""
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from exp46_dfg_build import LAD, NAMES, CACHE
from openfoam_cylinder import DFGCylinder, REF
from knobs import LocalBudget, run_protocol, screen
from replay import grid, method_cell, baseline_cell, fmt_method
from exp47_transfer import failure_rate, PHIS

KINDS = ["discretisation", "nested-solver", "nested-solver", "cost-only"]
HS = [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), [1.0 / c for c in LAD[3]]]
raw = pickle.load(open(CACHE, "rb"))
tensor = {i: (v[0], v[1], v[2]) for i, v in raw.items()}
ALL = grid(LAD)
TUT = lambda na: (LAD[0].index(na), 0, LAD[2].index(1e-6), 0)   # tutorial solver settings
show = lambda i: "n_arc=%d, res_ctrl=%.0e, lin_tol=%.0e, n_coarsest=%d" % tuple(
    LAD[k][j] for k, j in enumerate(i))
rel = lambda i, ob: abs(tensor[i][0][ob] - REF[ob]) / abs(REF[ob])

print("=" * 104)
print("OpenFOAM 1912 simpleFoam, DFG flow around a cylinder 2D-1 (Re = 20); "
      f"{len(ALL)} configurations")
print("reference values (Schaefer-Turek; John-Matthies):", REF)
nconv = sum(v[2].get("converged", False) for v in tensor.values())
print(f"runs meeting res_ctrl before the iteration limit: {nconv}/{len(ALL)}")
print("=" * 104)

fin = (len(LAD[1]) - 1, len(LAD[2]) - 1, 0)
print("\nrelative error at res_ctrl = 1e-7, lin_tol = 1e-8, against n_arc "
      "(cells: " + ", ".join(str(tensor[(j,) + fin][2].get("cells")) for j in range(len(LAD[0]))) + ")")
for ob in REF:
    e = [rel((j,) + fin, ob) for j in range(len(LAD[0]))]
    print(f"  {ob:<3} " + "  ".join(f"{x:.2e}" for x in e))
    print(f"  {'ord':<3} " + "  ".join(f"{np.log(e[j]/e[j+1])/np.log(LAD[0][j+1]/LAD[0][j]):8.2f}"
                                      for j in range(len(e) - 1)))

print("\nthe OpenFOAM tutorial solver settings (res_ctrl 1e-2/1e-3, GAMG 1e-6, "
      "coarsest 10) against tight ones, by mesh level:")
for na in LAD[0]:
    i = TUT(na)
    tight = (i[0], len(LAD[1]) - 1, len(LAD[2]) - 1, 0)
    print(f"  n_arc {na:>3}: " + "   ".join(
        f"{ob} {rel(i, ob):.2e} (solver part {abs(tensor[i][0][ob]-tensor[tight][0][ob])/abs(REF[ob]):.1e})"
        for ob in REF) + f"   cost {tensor[i][1]:.1f}s vs {tensor[tight][1]:.1f}s")

def get(idx, ob):
    v = tensor[idx]
    return v[0][ob] if np.isfinite(v[1]) else np.nan
get.exact = REF
print("\ntransfer of thresholds (Table 3 on OpenFOAM), pooled over cd and dp:")
print(f"  {'parameter':<12} " + " ".join(f"phi={p:<6g}" for p in PHIS))
for k in range(4):
    row = []
    for phi in PHIS:
        f, n = failure_rate(get, LAD, k, ["cd", "dp"], phi)
        row.append(f"{100*f/max(n,1):>9.0f}%")
    print(f"  {NAMES[k]:<12} " + " ".join(row))

TARGETS = {"cd": [1e-2, 3e-3, 1e-3], "dp": [3e-2, 1e-2], "cl": [1e-1, 3e-2]}
for ob, tgs in TARGETS.items():
    print(f"\nerror budget at the tutorial settings, n_arc = 16, for {ob}:")
    look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    lb = LocalBudget(look, LAD, TUT(16), hs=HS)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    tot = sum(abs(x) for x in lb.rk if np.isfinite(x)) or 1
    print("  " + ", ".join(f"{NAMES[k]} {100*abs(lb.rk[k])/tot:.0f}%" for k in range(4))
          + f";  estimate {abs(lb.E_knob)/abs(REF[ob]):.2e}, true {rel(TUT(16), ob):.2e}")
    for tg in tgs:
        eps = tg * abs(REF[ob])
        feas = [i for i in ALL if np.isfinite(tensor[i][1]) and rel(i, ob) <= tg]
        if not feas:
            print(f"  target {tg:.0e}: unattainable on the grid")
            continue
        print(f"\n--- {ob}, target {tg:.0e} relative")
        o = method_cell(tensor, LAD, HS, KINDS, ob, REF[ob], eps)
        print(f"  {'method':<28} " + fmt_method(o))
        for m in o["misses"][:6]:
            print(f"      miss: {show(m[0])} -> {show(m[1])}: true {m[2]:.2f} eps, estimate {m[3]:.2f} eps")
        o2 = method_cell(tensor, LAD, HS, KINDS, ob, REF[ob], eps, resid_margin=5.0)
        print(f"  {'method, gated margin':<28} " + fmt_method(o2))
        for lab, r in baseline_cell(tensor, LAD, ob, REF[ob], eps).items():
            print(f"  {lab:<28} attained {r['attained']}/{r['n']}, penalty median "
                  f"{r['pen']:.2f}, runs median {r['runs']:.0f}")
        best = o["opt"]
        for na in (16, 32):
            r = run_protocol(look, LAD, TUT(na), eps, hs=HS, kinds=KINDS)
            q = r["proposal"]
            print(f"  from the tutorial settings at n_arc={na}: {rel(TUT(na), ob):.2e} "
                  f"({tensor[TUT(na)][1]:.1f}s) -> {'certified' if r['certified'] else 'declined'} "
                  f"{show(q)}: {rel(q, ob):.2e} ({tensor[q][1]:.1f}s), {r['runs']} runs "
                  f"[cheapest attaining: {show(best)}, {tensor[best][1]:.1f}s]")

# a false declaration, measured live
print("\nrelax_U declared cost-only (it is documented as affecting only convergence):")
p = DFGCylinder()
L5 = LAD + [[0.9, 0.7, 0.5]]
store = {}
def live(idx):
    idx = tuple(idx)
    if idx not in store:
        th = {n: L5[k][j] for k, (n, j) in enumerate(zip(NAMES + ["relax_U"], idx))}
        obs, c, w, info = p.run(th)
        store[idx] = (obs["cd"], c)
    return store[idx]
anchor = (LAD[0].index(24), 4, 3, 0, 0)
active, frozen, viol, sc = screen(live, L5, anchor, KINDS + ["cost-only"])
print(f"  screen at {sc}: violations reported for "
      f"{[(NAMES + ['relax_U'])[k] for k in viol]}; frozen {[(NAMES + ['relax_U'])[k] for k in frozen]}")
for j, rx in enumerate(L5[4]):
    print(f"  relax_U = {rx}: c_D = {live(anchor[:4] + (j,))[0]:.6f}")
