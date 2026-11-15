"""Sensitivity of the round-off measurement of GROMACS (S-5.4).  The six settings of
exp44_gromacs_floor.py give six potential energies per build; for every subset of two to
six of them and every factor F in {3, 5, 10, 20}, the program declines targets below F
times their spread, and the mixed-precision runs are searched again.  The cost of the
three measurement runs is timed against runs at the GROMACS default setting; without the
GROMACS binaries the ratio measured before (figures/rev3_floor_sensitivity.txt) is used.

    python exp58_floor_sensitivity.py v1.2.1   # prints the table, writes figures/tab_floor.tex
"""
import itertools
import os
import sys

import numpy as np

import cases as C
import pooled as P

VALS = {"gmx": [-73672.5625, -73672.4297, -73672.8672, -73671.8906, -73670.6719, -73672.8672],
        "gmx_d": [-73671.3702, -73671.3702, -73671.3702, -73671.3702, -73671.3701, -73671.3702]}
CFG = sys.argv[1] if len(sys.argv) > 1 else "v1.2"
case = [c for c in C._gromacs() if c.key == "GMX-s"][0]
QD = case.exact["potential"]
v = np.array(VALS["gmx"])
print(f"configuration {CFG}; mixed-precision values of the six converged settings: {v.tolist()}")
TARGETS = [1e-4, 3e-5, 1e-5]


def sci(x):
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(x) + 1e-12))
    m = x / 10**e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.1f}{{\\times}}10^{{{e}}}$"


TEX = ["\\begin{tabular}{rrrlrrrr}", "\\toprule",
       "& & & & \\multicolumn{2}{c}{Accepted at $10^{-4}$} & Accepted & \\\\", "\\cmidrule(lr){5-6}",
       "$k$ & $F$ & Floors & Declines below (relative) & smallest & largest & tighter & Missed \\\\",
       "\\midrule"]
for k in range(2, 7):
    floors = sorted({float(np.ptp(v[list(s)])) for s in itertools.combinations(range(6), k)})
    for F in (3, 5, 10, 20):
        thr = [F * f / abs(QD) for f in floors]
        line = (f"k={k} F={F:<2}: declines below {min(thr):.1e}..{max(thr):.1e} relative "
                f"({len(floors)} distinct floors)")
        cells = []
        allres = []
        for tg in TARGETS:
            res = []
            for f in (min(floors), max(floors)):
                case.floor = f
                cfg = dict(P.CONFIGS[CFG])
                kw = dict(cfg["kw"]); kw["floor_factor"] = F
                cfg = dict(cfg, kw=kw)
                o = P.replay(case, "potential", tg, cfg)
                res.append((o["accepted"], o["missed"]))
            cells.append(f"{tg:.0e}: acc {res[0][0]}-{res[1][0]}, miss {res[0][1]}-{res[1][1]}")
            allres.append(res)
        print(line + " | " + "; ".join(cells), flush=True)
        lo, hi = min(thr), max(thr)
        rng = sci(hi) if lo == hi else f"{sci(lo)}--{sci(hi)}"
        tighter = sum(r[0] for res in allres[1:] for r in res)
        missed = sum(r[1] for res in allres for r in res)
        TEX.append(f"{k} & {F} & {len(floors)} & {rng} & {allres[0][0][0]} & {allres[0][1][0]} & "
                   f"{tighter} & {missed} \\\\")
TEX += ["\\bottomrule", "\\end{tabular}"]
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"), exist_ok=True)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures", "tab_floor.tex"), "w") as fh:
    fh.write("\n".join(TEX) + "\n")

# Cost of the measurement.  The three measurement settings and the GROMACS
# default setting are timed together now, so that machine load cancels in their
# ratio; the ratio is applied to the tensor's cost of the default setting and set
# against the calibration cost of the searches (sum of the costs of the runs each
# search requests), median over all starts.
import re
import shutil
from gromacs_water import GromacsWater
from errbudget import run_protocol

STORED = {}
_rev3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures", "rev3_floor_sensitivity.txt")
if os.path.exists(_rev3):
    for _line in open(_rev3):
        _m = re.match(r"(gmx_d|gmx): three probe runs cost ([0-9.]+) default runs", _line)
        if _m:
            STORED[_m.group(1)] = float(_m.group(2))
CONVERGED = [dict(rcoulomb=rc, fourierspacing=fs, pme_order=po, ewald_rtol=rt)   # as exp44
             for rc, fs, po, rt in [(1.4, 0.06, 8, 1e-8), (1.4, 0.05, 8, 1e-9), (1.5, 0.06, 8, 1e-9)]]
DEFAULT = dict(rcoulomb=1.0, fourierspacing=0.12, pme_order=4, ewald_rtol=1e-5)
didx = (1, 2, 0, 1)
for build, key in (("gmx_d", "GMX-d"), ("gmx", "GMX-s")):
    cs = [c for c in C._gromacs() if c.key == key][0]
    if shutil.which(build):
        g = GromacsWater(nbox=2, nsteps=40, binary=build)
        td = [g.run(DEFAULT)[1] for _ in range(3)]
        tp = [g.run(th)[1] for th in CONVERGED[:3]]
        ratio = sum(tp) / np.median(td)
    elif build in STORED:
        ratio = STORED[build]
        print(f"{build}: binary not found; using the timing ratio of revision 3 "
              f"({ratio:.1f} default runs for the three probe runs)", flush=True)
    else:
        print(f"{build}: binary not found and no stored timing ratio; probe cost skipped", flush=True)
        continue
    probe_cost = ratio * cs.tensor[didx][1]
    tg = cs.targets["potential"][0]
    eps = tg * abs(cs.exact["potential"])
    look = lambda i: ((cs.tensor[tuple(i)][0]["potential"], cs.tensor[tuple(i)][1])
                      if np.isfinite(cs.tensor[tuple(i)][1]) else (np.nan, np.inf))
    cal = []
    for a in cs.grid():
        if not cs.usable(a):
            continue
        r = run_protocol(look, cs.LAD, a, eps, hs=cs.HS, kinds=cs.KINDS, **P.CONFIGS[CFG]["kw"])
        cal.append(r["protocol_cost"])
    print(f"{build}: three probe runs cost {ratio:.1f} default runs ({probe_cost*1e3:.1f} ms/step in all); "
          f"calibration of one search at {tg:.0e}: median {np.median(cal)*1e3:.1f} ms/step "
          f"({np.median(cal)/cs.tensor[didx][1]:.0f} default runs); probe/calibration "
          f"{probe_cost/np.median(cal):.2f}", flush=True)
