"""Item: what does each wrong class declaration cost?  (L3, all 432 starts)

Each of the five L3 parameters is declared as each of the two classes it does
not belong to (10 misdeclarations), and the method is replayed from every
start at 1e-2 and 3e-3 relative for both observables.  Declarations act in
the screen (which parameters may be frozen, and at which setting influence is
measured) and in the cost-only pass, so the screened mode -- the default of
knobs.tune -- is the one that can be affected.  The full mode is reported too.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import grids
from replay import method_cell

P = grids.L3
t = grids.load(P)
EX = grids.exact(P)
HS = grids.hs(P)
TRUE = P["kinds"]
CLASSES = ["discretisation", "nested-solver", "cost-only"]
CONF = [("correct", TRUE)]
for k, nm in enumerate(P["names"]):
    for c in CLASSES:
        if c != TRUE[k]:
            kk = list(TRUE); kk[k] = c
            CONF.append((f"{nm}: {TRUE[k]} -> {c}", kk))

for mode, scr in [("screened (default)", True), ("full", False)]:
    print("=" * 110)
    print(f"L3, {mode} mode: pooled over energy/gradient x 1e-2/3e-3, 432 starts each")
    print(f"  {'declaration':<42} {'certified':>10} {'missed':>7} {'safety':>8} "
          f"{'med pen':>8} {'worst pen':>9} {'med runs':>9}")
    for lab, kinds in CONF:
        cert = miss = 0; pen = []; runs = []
        for ob in P["obs"]:
            for tg in P["targets"]:
                o = method_cell(t, P["lad"], HS, kinds, ob, EX[ob], tg * abs(EX[ob]),
                                use_screen=scr)
                cert += o["certified"]; miss += o["missed"]
                pen += o["pen"]; runs += o["runs"]
        print(f"  {lab:<42} {cert:>10} {miss:>7} {100*(cert-miss)/max(cert,1):>7.1f}% "
              f"{np.median(pen):>8.2f} {np.max(pen):>9.2f} {np.median(runs):>9.0f}")

# ---------------------------------------------------------------------------
# The same test on OpenFOAM (DFG 2D-1), where a nested-solver tolerance does
# carry a large share of the error (the lift).  Screened mode only.
import pickle as _pk
from exp46_dfg_build import LAD as DLAD, NAMES as DN, CACHE as DC
from openfoam_cylinder import REF as DREF
if os.path.exists(DC):
    dt = _pk.load(open(DC, "rb"))
    DT = ["discretisation", "nested-solver", "nested-solver", "cost-only"]
    DHS = [[1.0 / n for n in DLAD[0]], list(DLAD[1]), list(DLAD[2]), [1.0 / c for c in DLAD[3]]]
    conf = [("correct", DT)]
    for k, nm in enumerate(DN):
        for c in CLASSES:
            if c != DT[k]:
                kk = list(DT); kk[k] = c
                conf.append((f"{nm}: {DT[k]} -> {c}", kk))
    print("=" * 110)
    print("OpenFOAM DFG 2D-1, screened mode: c_D at 3e-3 and c_L at 3e-2, 780 starts each")
    print(f"  {'declaration':<42} {'certified':>10} {'missed':>7} {'safety':>8} "
          f"{'med pen':>8} {'worst pen':>9} {'med runs':>9}")
    for lab, kinds in conf:
        cert = miss = 0; pen = []; runs = []
        for ob, tg in [("cd", 3e-3), ("cl", 3e-2)]:
            o = method_cell(dt, DLAD, DHS, kinds, ob, DREF[ob], tg * abs(DREF[ob]),
                            use_screen=True)
            cert += o["certified"]; miss += o["missed"]; pen += o["pen"]; runs += o["runs"]
        print(f"  {lab:<42} {cert:>10} {miss:>7} {100*(cert-miss)/max(cert,1):>7.1f}% "
              f"{np.median(pen):>8.2f} {np.max(pen):>9.2f} {np.median(runs):>9.0f}")
