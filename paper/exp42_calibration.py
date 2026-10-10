"""Section 7.7 and Table 10: calibration cost, saving, break-even, screening,
and calibration runs against the ladder distance between start and answer (L3)."""
import sys, os, collections
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from knobs import run_protocol
from grids import L3, hs, load, exact, grid

t, ex, H, ALL = load(L3), exact(L3), hs(L3), grid(L3)
print(f"{'observable':<10} {'target':<7} {'mode':<9} {'K':>7} {'runs':>5} {'s':>6} {'n_BE':>6} {'penalty':>8}")
res = {}
for ob in L3["obs"]:
    look = lambda i, ob=ob: (t[tuple(i)][0][ob], t[tuple(i)][1])
    for tg in L3["targets"]:
        eps = tg * abs(ex[ob])
        feas = [i for i in ALL if abs(t[i][0][ob] - ex[ob]) <= eps]
        copt = min(t[i][1] for i in feas)
        for use in (False, True):
            K, runs, sav, pen = [], [], [], []
            for a in ALL:
                r = run_protocol(look, L3["lad"], a, eps, hs=H, kinds=L3["kinds"], use_screen=use)
                q = r["proposal"]
                K.append(r["protocol_cost"] / t[a][1]); runs.append(r["runs"])
                pen.append(t[q][1] / copt)
                if abs(t[a][0][ob] - ex[ob]) <= eps:
                    sav.append(t[a][1] / t[q][1])
            K, s = float(np.median(K)), float(np.median(sav))
            mode = "screened" if use else "full"
            res[(ob, tg, mode)] = (K, s)
            print(f"{ob:<10} {tg:<7.0e} {mode:<9} {K:>7.1f} {np.median(runs):>5.0f} {s:>6.2f}"
                  f" {K/(1-1/s):>6.0f} {np.median(pen):>8.2f}")
print("\nproduction runs below which screening is cheaper overall "
      "(calibration plus production, in units of the user's run):")
for ob in L3["obs"]:
    for tg in L3["targets"]:
        Kf, sf = res[(ob, tg, "full")]; Ks, ss = res[(ob, tg, "screened")]
        print(f"  {ob:<10} {tg:.0e}: {(Kf-Ks)/(1/ss-1/sf):.0f}")

print("\ncalibration runs against the ladder distance |start - answer|_1 (L3, all cells):")
by = collections.defaultdict(list)
for ob in L3["obs"]:
    look = lambda i, ob=ob: (t[tuple(i)][0][ob], t[tuple(i)][1])
    for tg in L3["targets"]:
        eps = tg * abs(ex[ob])
        for a in ALL:
            r = run_protocol(look, L3["lad"], a, eps, hs=H, kinds=L3["kinds"], use_screen=False)
            by[sum(abs(x - y) for x, y in zip(a, r["proposal"]))].append(r["runs"])
print("  " + ", ".join(f"{d}: {int(np.median(v))}" for d, v in sorted(by.items())))
