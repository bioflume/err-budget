"""Where is the cliff in the safety margin, and how wide is the safe plateau?

The referee's objection is exact: F_s = 1.25 was chosen on the same cells on
which reliability is then reported, so the reported reliability restates the
selection criterion.  This answers it in the only way that is worth anything:

  * calibrate on L2 alone and hold L3 and L4 out;
  * scan F_s finely enough to locate the transition rather than bracket it;
  * report, per cell, the SMALLEST F_s on the scan that was safe, so a reader
    can see the headroom instead of taking 1.25 on trust;
  * report the cost of the margin at the same time, because a margin that buys
    safety by refining everything is not a result.
"""
import itertools
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget import run_protocol
from errbudget.problems import Burgers1D, Burgers2D, ReactionDiffusion2D

D = os.path.join(os.path.dirname(__file__), "data")
FS = [1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.35, 1.50, 2.00, 3.00]


def scan(tensor, LAD, HS, KIN, ob, QS, t, rel):
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    eps = t * abs(QS) if rel else t
    err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
    feas = [i for i in ALL if err[i] <= eps]
    if not feas:
        return None
    c_opt = min(tensor[i][1] for i in feas)
    lk = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    row = {}
    for fs in FS:
        safe = cert = 0
        pens = []
        for a in ALL:
            r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KIN,
                             use_screen=False, fs_cap=fs)
            q = r["proposal"]
            if r["certified"]:
                cert += 1
                safe += err[q] <= eps
                pens.append(tensor[q][1] / c_opt)
        row[fs] = dict(certified=cert, safe=safe, n=len(ALL),
                       pen_med=float(np.median(pens)) if pens else np.nan,
                       pen_max=float(np.max(pens)) if pens else np.nan)
    return row


b1 = Burgers1D(nu=0.05, T=1.0)
t1 = pickle.load(open(os.path.join(D, "burgers_tensor.pkl"), "rb"))
LAD1 = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
        [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
HS1 = [[2 * np.pi / n for n in LAD1[0]], [1.0 / s for s in LAD1[1]],
       list(LAD1[2]), list(LAD1[3])]
KIN1 = ["discretisation", "discretisation", "nested-solver", "nested-solver"]

b2 = Burgers2D(nu=0.05, T=1.0)
t2 = pickle.load(open(os.path.join(D, "burgers2d_tensor.pkl"), "rb"))
LAD2 = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS2 = [[2 * np.pi / n for n in LAD2[0]], [1.0 / s for s in LAD2[1]],
       list(LAD2[2]), list(LAD2[3]), list(LAD2[4])]
KIN2 = ["discretisation", "discretisation", "nested-solver", "nested-solver",
        "cost-only"]

p4 = ReactionDiffusion2D()
t4 = pickle.load(open(os.path.join(D, "reactdiff_tensor.pkl"), "rb"))
LAD4 = [[16, 32, 64, 128], [5, 10, 20, 40], [1, 2, 4, 8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS4 = [[2 * np.pi / n for n in LAD4[0]], [1.0 / s for s in LAD4[1]],
       [1.0 / r for r in LAD4[2]], list(LAD4[3]), list(LAD4[4])]
KIN4 = ["discretisation", "discretisation", "discretisation",
        "nested-solver", "cost-only"]

CELLS = [
    ("L2", "calibration", t1, LAD1, HS1, KIN1, b1.exact_observables(),
     ["energy", "dissipation"], [1e-2, 3e-3, 1e-3], False),
    ("L3", "held out", t2, LAD2, HS2, KIN2, b2.exact_observables(),
     ["energy", "gradient"], [1e-2, 3e-3], True),
    ("L4", "held out", t4, LAD4, HS4, KIN4, p4.exact_observables(),
     ["inhibitor_mass", "pattern_sharpness"], [1e-2, 3e-3], True),
]

results = []
print("=" * 108)
print("safety and cost against the margin F_s; F_s is fixed on L2 and L3/L4 "
      "are held out")
print("=" * 108)
hdr = "  ".join(f"{f:>5.2f}" for f in FS)
print(f"{'cell':<34} {'role':<12} {hdr}")
for name, role, T, LAD, HS, KIN, EX, obs, targets, rel in CELLS:
    for ob in obs:
        for t in targets:
            row = scan(T, LAD, HS, KIN, ob, EX[ob], t, rel)
            if row is None:
                continue
            lab = f"{name} {ob} {t:.0e}" + ("rel" if rel else "")
            safe_txt = "  ".join(
                f"{100*row[f]['safe']/max(row[f]['certified'],1):>5.0f}"
                for f in FS)
            print(f"{lab:<34} {role:<12} {safe_txt}   % safe")
            pen_txt = "  ".join(f"{row[f]['pen_max']:>5.1f}" for f in FS)
            print(f"{'':<34} {'':<12} {pen_txt}   worst-case cost penalty")
            ok = [f for f in FS
                  if row[f]["safe"] == row[f]["certified"] == row[f]["n"]]
            fmin = min(ok) if ok else None
            print(f"{'':<34} {'':<12} smallest fully safe margin on this scan: "
                  f"{fmin if fmin else 'none'}")
            results.append(dict(cell=lab, role=role, fmin=fmin,
                                row={str(f): row[f] for f in FS}))

json.dump(results, open(os.path.join(D, "margin_scan.json"), "w"), indent=1)

print()
print("=" * 108)
cal = [r for r in results if r["role"] == "calibration"]
held = [r for r in results if r["role"] == "held out"]
fc = [r["fmin"] for r in cal if r["fmin"]]
fh = [r["fmin"] for r in held if r["fmin"]]
print(f"  calibration cells: smallest fully safe margin ranges "
      f"{min(fc):.2f}-{max(fc):.2f}; choosing the largest gives F_s = {max(fc):.2f}")
print(f"  held-out cells   : smallest fully safe margin ranges "
      f"{min(fh):.2f}-{max(fh):.2f}")
print(f"  the margin chosen on L2 is {'SUFFICIENT' if max(fh) <= max(fc) else 'INSUFFICIENT'}"
      f" for every held-out cell")
print(f"  headroom at F_s = 1.25: the worst held-out cell needed {max(fh):.2f}")
