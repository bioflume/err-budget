"""Sensitivity of the GROMACS round-off floor to the number of probe settings and
to the factor (referee points 22 and 23 of revision 3).

The six converged settings of exp44_gromacs_floor.py give six values of the
potential energy per build (figures/sec82_floor.txt).  For every subset of k of
them (k = 2..6) the floor is their spread; for each factor F the program
declines targets below F times the floor.  For each (k, F) we report, over all
subsets, the range of the declining threshold and, replaying the mixed-precision
tensor with that floor, the accepted starts and the misses against the
double-precision reference.  The probe's cost is reported relative to the
calibration cost of one search.
"""
import itertools
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C
import pooled as P

VALS = {"gmx": [-73672.5625, -73672.4297, -73672.8672, -73671.8906, -73670.6719, -73672.8672],
        "gmx_d": [-73671.3702, -73671.3702, -73671.3702, -73671.3702, -73671.3701, -73671.3702]}
CFG = sys.argv[1] if len(sys.argv) > 1 else "v1.2"
case = [c for c in C.all_cases() if c.key == "GMX-s"][0]
QD = case.exact["potential"]
v = np.array(VALS["gmx"])
print(f"configuration {CFG}; mixed-precision values of the six converged settings: {v.tolist()}")
TARGETS = [1e-4, 3e-5, 1e-5]
for k in range(2, 7):
    floors = sorted({float(np.ptp(v[list(s)])) for s in itertools.combinations(range(6), k)})
    for F in (3, 5, 10, 20):
        thr = [F * f / abs(QD) for f in floors]
        line = (f"k={k} F={F:<2}: declines below {min(thr):.1e}..{max(thr):.1e} relative "
                f"({len(floors)} distinct floors)")
        cells = []
        for tg in TARGETS:
            acc_lo = miss_hi = None
            res = []
            for f in (min(floors), max(floors)):
                case.floor = f
                cfg = dict(P.CONFIGS[CFG])
                kw = dict(cfg["kw"]); kw["floor_factor"] = F
                cfg = dict(cfg, kw=kw)
                o = P.replay(case, "potential", tg, cfg)
                res.append((o["accepted"], o["missed"]))
            cells.append(f"{tg:.0e}: acc {res[0][0]}-{res[1][0]}, miss {res[0][1]}-{res[1][1]}")
        print(line + " | " + "; ".join(cells), flush=True)

# cost of the probe: six converged settings against the calibration of one search
t, _ = pickle.load(open(os.path.join(C.D, "gromacs_tensor_gmx.pkl"), "rb"))
cmax = max(x[1] for x in t.values() if np.isfinite(x[1]))
print(f"\nthe most expensive setting on the GROMACS grid costs {cmax*1000:.1f} ms/step; "
      "the probe settings (rcoulomb 1.4-1.6 nm, spacing 0.05-0.06 nm, order 8-10) cost "
      "about 2-4 times that each (exp44 timings)")
