"""Quantum ESPRESSO forces and energies against conv_thr (item 6).

Exhaustive tensor: ecutwfc 8 x nk 4 x conv_thr 5 x mixing_beta 3 = 480 runs,
graded against a reference beyond the ladder in every parameter.  Reports how
much of each observable's error budget the SCF threshold carries, whether a
threshold adequate for the energy is adequate for the force, and the
method's replay from every start (power-law and exponential cutoff models).
"""
import os
import pickle
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from qe_silicon import QESilicon
from knobs import LocalBudget
from replay import grid, method_cell, baseline_cell, fmt_method

NAMES = ["ecutwfc", "nk", "conv_thr", "mixing_beta"]
LAD = [[12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 40.0, 48.0], [2, 3, 4, 6],
       [1e-4, 1e-6, 1e-8, 1e-10, 1e-12], [0.3, 0.5, 0.7]]
KINDS = ["discretisation", "discretisation", "nested-solver", "cost-only"]
HS = [[1 / np.sqrt(e) for e in LAD[0]], [1.0 / k for k in LAD[1]], list(LAD[2]),
      [1.0 / b for b in LAD[3]]]
HS_EXP = [list(LAD[0])] + HS[1:]
# The exponential model extrapolates further than the power law; its own fit
# residual sets the margin (the residual-gated margin of Section 7.4, alpha = 5).
MODELS = {"power": (HS, None, {}),
          "exponential": (HS_EXP, ["exponential", "power", "power", "power"], {}),
          "exp. + gated": (HS_EXP, ["exponential", "power", "power", "power"],
                           dict(resid_margin=5.0))}
CACHE = os.path.join(os.path.dirname(__file__), "data", "qe_tensor.pkl")

p = QESilicon()
if os.path.exists(CACHE):
    tensor, REF, UNC = pickle.load(open(CACHE, "rb"))
else:
    t0 = time.time()
    part = CACHE + ".partial"
    saved = pickle.load(open(part, "rb")) if os.path.exists(part) else {}
    def once(key, th):
        if key not in saved:
            saved[key] = p.run(th)[0]
            pickle.dump(saved, open(part, "wb"))
        return saved[key]
    REF = once("ref", p.reference())
    # uncertainty: last cutoff increment at the reference k-mesh, plus the
    # last k-mesh increment at a lower cutoff (the two are nearly separable)
    e75 = once("e75", dict(p.reference(), ecutwfc=75.0))
    k8 = once("k8", dict(p.reference(), ecutwfc=60.0))
    k10 = once("k10", dict(p.reference(), ecutwfc=60.0, nk=10))
    UNC = {ob: (abs(e75[ob] - REF[ob]) + abs(k10[ob] - k8[ob])) / abs(REF[ob]) for ob in REF}
    print(f"reference {REF}; uncertainty {UNC} ({time.time()-t0:.0f}s)", flush=True)
    tensor = saved.setdefault("tensor", {})
    for n, idx in enumerate(grid(LAD)):
        if idx in tensor:
            continue
        th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
        obs, c, w, info = p.run(th)
        tensor[idx] = (obs, c, info)
        if n % 10 == 0:
            pickle.dump(saved, open(part, "wb"))
            print(f"  {n}/{len(grid(LAD))}  {time.time()-t0:.0f}s", flush=True)
    pickle.dump((tensor, REF, UNC), open(CACHE, "wb"))
    os.remove(part)
    print(f"tensor: {len(tensor)} pw.x runs in {time.time()-t0:.0f}s")

ok = [i for i in grid(LAD) if np.isfinite(tensor[i][1])]
print("=" * 100)
print(f"Quantum ESPRESSO 6.7, Si with a displaced atom; {len(ok)}/{len(tensor)} runs usable")
print(f"reference: E = {REF['energy']:.10f} Ry, |F| = {REF['force']:.8e} Ry/bohr; "
      f"relative uncertainty {UNC['energy']:.1e} (E), {UNC['force']:.1e} (F)")
print("=" * 100)

print("\nrelative error against conv_thr at ecutwfc = 48 Ry, 6x6x6 k-points, beta = 0.7:")
for ob in ["energy", "force"]:
    e = [abs(tensor[(7, 3, j, 2)][0][ob] - REF[ob]) / abs(REF[ob]) for j in range(5)]
    d = [abs(tensor[(7, 3, j, 2)][0][ob] - tensor[(7, 3, 4, 2)][0][ob]) / abs(REF[ob]) for j in range(5)]
    print(f"  {ob:<7} total " + "  ".join(f"{x:.1e}" for x in e))
    print(f"  {'':<7} from conv_thr alone " + "  ".join(f"{x:.1e}" for x in d))

print("\nloosest conv_thr whose own contribution stays below 10% of the total error, "
      "by ecutwfc and k-points (beta = 0.7):")
for ob in ["energy", "force"]:
    print(f"  {ob}:")
    for i in range(len(LAD[0])):
        row = []
        for k in range(len(LAD[1])):
            fine = tensor[(i, k, 4, 2)][0][ob]
            tot = abs(fine - REF[ob])
            j = next(j for j in range(5) if abs(tensor[(i, k, j, 2)][0][ob] - fine) <= 0.1 * tot)
            row.append(f"{LAD[2][j]:.0e}")
        print(f"    ecut {LAD[0][i]:4.0f} Ry: " + "  ".join(row) + "   (nk = 2, 3, 4, 6)")

TARGETS = {"energy": [1e-3, 1e-4, 1e-5], "force": [1e-1, 3e-2, 1e-2]}
START = (1, 1, 1, 2)    # ecut 16 Ry, 3x3x3, conv_thr 1e-6 (the pw.x default), beta 0.7
for ob, tgs in TARGETS.items():
    for tg in tgs:
        if tg < 10 * UNC[ob]:
            print(f"\n--- {ob} {tg:.0e}: below ten times the reference uncertainty; skipped")
            continue
        eps = tg * abs(REF[ob])
        print(f"\n--- {ob}, target {tg:.0e} relative")
        for mname, (hs, models, kw) in MODELS.items():
            o = method_cell(tensor, LAD, hs, KINDS, ob, REF[ob], eps, models=models, **kw)
            print(f"  method, {mname:<12} " + fmt_method(o))
        for lab, r in baseline_cell(tensor, LAD, ob, REF[ob], eps).items():
            print(f"  {lab:<28} attained {r['attained']}/{r['n']}, penalty median "
                  f"{r['pen']:.2f}, runs median {r['runs']:.0f}")
        look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
        lb = LocalBudget(look, LAD, START, hs=HS_EXP, models=MODELS["exponential"][1])
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        tot = sum(abs(x) for x in lb.rk if np.isfinite(x)) or 1
        print("  budget at the pw.x-default start: " + ", ".join(
            f"{NAMES[k]} {100*abs(lb.rk[k])/tot:.0f}%" for k in range(4)))
