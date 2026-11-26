"""Numbers of the supplement that no script printed (REPRODUCTION_CHECK.md, "Not reproduced").

Each part recomputes one group of numbers from the stored runs and prints it next to the
value the supplement quotes, so that a reader can see whether it is reproduced.

1. S-5.1: the shares of momentum (U) solves stopped by their tolerance alone, over the 608
   converged runs of the 2D-1 grid with OpenFOAM v2406 (26% at 100 t/r_c = 0.1, 94% at 1).
   "Tolerance alone" is the classes abs (stopped on the absolute tolerance) and zero (no
   iteration: the initial residual was already below it) of
   ../hamilton/J1_dfg2d1/analyse_j1.py, which prints the same shares over all 840 runs.
2. S-5.2: the largest Courant number of the refined 2D-2 runs, Co_max = 27.8 n_arc dt, and
   "22 of 54 runs diverged".
3. S-3.2: the LAMMPS errors along the five cutoffs, (3.1, 3.6, 0.77, 2.4, 0.95)e-2, and the
   ABINIT estimate at the finest setting (1.15 and 3.5 times the target).

    python exp85_unscripted.py      # seconds; prints only, as in figures/unscripted.txt
"""
import os
import pickle
from collections import Counter

import numpy as np

import cases as C
import pooled as P
from exp46_dfg_build import LAD as LAD_2D1
from paths import HERE

RESULTS = os.path.join(HERE, "..", "hamilton", "results")


def momentum_stops():
    d = pickle.load(open(os.path.join(RESULTS, "dfg_tensor_current.pkl"), "rb"))
    print("1. S-5.1: momentum solves stopped by their tolerance alone (abs + zero), last 20"
          " outer iterations")
    print(f"   {'100 t/r_c':>9} {'runs':>12} {'all runs':>9} {'converged':>10}  supplement")
    quoted = {-1: "26%", 0: "94%"}
    for key in (-1, 0):
        out = {}
        for only_converged in (False, True):
            counts, runs = Counter(), 0
            for i, (_, _, info) in d.items():
                stop = info.get("stop")
                if not stop or (only_converged and not info.get("converged")):
                    continue
                rc, t = LAD_2D1[1][i[1]], LAD_2D1[2][i[2]]
                # U is solved to 10 t against residualControl rc / 10 (analyse_j1.stop_summary)
                if int(round(np.log10(100 * t / rc))) != key:
                    continue
                runs += 1
                counts.update(stop["last"]["U"])
            out[only_converged] = (runs, (counts["abs"] + counts["zero"]) / sum(counts.values()))
        n_conv = sum(1 for v in d.values() if v[2].get("converged"))
        print(f"   {10.0 ** key:>9g} {out[False][0]:>5} / {out[True][0]:<4}"
              f" {out[False][1]:9.1%} {out[True][1]:10.1%}  {quoted[key]}")
    print(f"   converged runs: {n_conv} of {len(d)}")


def courant():
    co = pickle.load(open(os.path.join(RESULTS, "j2_tensor.pkl"), "rb"))["courant"]
    rows = [(int(k.split("|")[0]), float(k.split("|")[1]), c) for k, v in co.items() for c in v]
    n, dt, c = (np.array(x) for x in zip(*rows))
    # a diverged run reports an overflowing Courant number (1e40 and more in j2.txt)
    diverged = c > 1e3
    ratio = c[~diverged] / (n[~diverged] * dt[~diverged])
    print("\n2. S-5.2: largest Courant number of the 2D-2 runs")
    print(f"   Co_max / (n_arc dt) over the {int((~diverged).sum())} runs that did not diverge:"
          f" median {np.median(ratio):.2f}, range {ratio.min():.2f}-{ratio.max():.2f}"
          f"   (supplement: 27.8)")
    print(f"   diverged: {int(diverged.sum())} of {len(c)} runs")
    partial = [(k, sum(x > 1e3 for x in v), len(v)) for k, v in co.items()
               if 0 < sum(x > 1e3 for x in v) < len(v)]
    print(f"   in the {len(partial)} (n_arc, dt) pairs where some but not all runs diverged: "
          f"{sum(p[1] for p in partial)} of {sum(p[2] for p in partial)} runs diverged"
          f"   (supplement: 22 of 54)")
    for k, nd, nt in sorted(partial, key=lambda p: (int(p[0].split("|")[0]), p[0])):
        print(f"     n_arc {k.split('|')[0]:>3}, dt {k.split('|')[1]:>7}: {nd} of {nt}")


def lammps_cutoffs():
    case = next(c for c in C._lammps() if c.key == "LMP-20")
    ob = "energy_per_ion"
    acc, order = 0, 3          # PPPM accuracy 1e-2, order 7
    errs = [abs(case.tensor[(k, acc, order)][0][ob] - case.exact[ob]) / abs(case.exact[ob])
            for k in range(len(case.LAD[0]))]
    print("\n3a. S-3.2: LAMMPS (8,000 ions), relative energy error along the cutoffs"
          f" {case.LAD[0]} at PPPM accuracy {case.LAD[1][acc]:g}, order {case.LAD[2][order]}")
    print("   " + ", ".join(f"{e * 100:.2g}" for e in errs) + " (x 1e-2)"
          "   (supplement: 3.1, 3.6, 0.77, 2.4, 0.95)")


def abinit_finest():
    case = C._abinit()
    finest = tuple(len(l) - 1 for l in case.LAD)
    kw = P.CONFIGS["v1.2.1"]["kw"]
    T = case.tensor
    print(f"\n3b. S-3.2: ABINIT, the program's error estimate at the finest setting {finest}"
          f" (v1.2.1), as a multiple of the target")
    print(f"   {'target':>7} {'estimate':>9} {'true error':>11} {'estimate/true':>14}  result")
    for tg in case.targets["energy"]:
        eps = tg * abs(case.exact["energy"])

        def look(i):
            v = T[tuple(i)]
            return (v[0]["energy"], v[1]) if np.isfinite(v[1]) else (np.nan, np.inf)

        r = P.run_protocol(look, case.LAD, finest, eps, hs=case.HS, kinds=case.KINDS,
                           models=case.models, **kw)
        est = abs(r["budget"].E0) / eps
        true = abs(T[finest][0]["energy"] - case.exact["energy"]) / eps
        res = "accepted" if r["certified"] else f"declined ({r.get('declined') or 'not_reached'})"
        print(f"   {tg:>7g} {est:9.3g} {true:11.3g} {est / true:14.3g}  {res}")
    print("   (supplement: 1.15 and 3.5 times the target; 'about 50 times')")


if __name__ == "__main__":
    momentum_stops()
    courant()
    lammps_cutoffs()
    abinit_finest()
