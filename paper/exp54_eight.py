"""Eight parameters on L3 (S-3.5): the five of L3 and gmres_restart, precond_freq and
newton_max, 11,664 settings that are not built.  The released program searches from 200
random starting settings at the energy target 1e-2, and every returned setting is graded
against the manufactured solution.  The stored runs (data/eight_cache.pkl) have no
convergence flags, so every setting a search requests is run once to obtain its flag and
checked against the stored run; the flags are cached in data/eight_flags_released.pkl.

    ERRBUDGET_NPROC=4 python exp54_eight.py   # 6.5 CPU-hours without the flag cache;
                                              # writes data/eight_released.json
"""
import os
import pickle
import time

import numpy as np

import json
from multiprocessing import Pool
from errbudget import run_protocol, V121
from paths import DATA, OUT
from errbudget.problems import Burgers2D
import grids

NAMES = ["N", "steps", "newton_tol", "gmres_tol", "ilu_drop", "gmres_restart",
         "precond_freq", "newton_max"]
LAD = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8], [1e-2, 1e-5, 1e-8],
       [1e-1, 1e-2, 1e-3], [10, 30, 100], [20, 5, 1], [2, 5, 20]]
KINDS = ["discretisation", "discretisation", "nested-solver", "nested-solver",
         "cost-only", "nested-solver", "cost-only", "nested-solver"]
HS = [[2 * np.pi / n for n in LAD[0]], [1.0 / s for s in LAD[1]], list(LAD[2]),
      list(LAD[3]), list(LAD[4]), [1.0 / r for r in LAD[5]], list(LAD[6]),
      [1.0 / m for m in LAD[7]]]
OB, TG = "energy", 1e-2
NSTART = int(os.environ.get("ERRBUDGET_NSTART", 200))
CACHE = os.path.join(DATA, "eight_cache.pkl")
FLAGS = os.path.join(OUT, "eight_flags_released.pkl")
NPROC = int(os.environ.get("ERRBUDGET_NPROC", 3))

b = Burgers2D(nu=0.05, T=1.0)
EX = b.exact_observables()
eps = TG * abs(EX[OB])
store = pickle.load(open(CACHE, "rb")) if os.path.exists(CACHE) else {}
flag = pickle.load(open(FLAGS, "rb")) if os.path.exists(FLAGS) else {}   # idx -> (obs, c, info)
t0, new = time.time(), [0]


def _run(idx):
    th = {NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}
    obs, c, w, info = Burgers2D(nu=0.05, T=1.0).run(th)
    return idx, (obs, c, dict(info))


def _record(idx, val):
    flag[idx] = val
    if idx in store:
        flag[idx][2]["same"] = (val[0] == store[idx][0] and val[1] == store[idx][1])
    else:
        store[idx] = (val[0], val[1])
    new[0] += 1
    if new[0] % 50 == 0:
        pickle.dump(flag, open(FLAGS, "wb"))
        print(f"    {len(flag)} flagged runs, {time.time()-t0:.0f}s", flush=True)


touched = set()


def look_dry(idx):
    """First pass: assume every run converged, record the settings requested."""
    idx = tuple(idx)
    touched.add(idx)
    if idx not in store:
        _record(*_run(idx))
    return store[idx][0][OB], store[idx][1]


def look(idx):
    idx = tuple(idx)
    if idx not in flag:
        _record(*_run(idx))
    obs, c, info = flag[idx]
    return (obs[OB] if info.get("converged", True) else np.nan), c


rng = np.random.default_rng(2026)
starts = sorted({tuple(int(rng.integers(0, len(l))) for l in LAD) for _ in range(3 * NSTART)})
starts = [starts[i] for i in rng.permutation(len(starts))[:NSTART]]
P5 = grids.L3
t5 = grids.load(P5)
feas5 = [i for i in grids.grid(P5) if abs(t5[i][0][OB] - EX[OB]) <= eps]
c5 = min(t5[i][1] for i in feas5)

def main():

    # first pass: which settings do the searches request?
    for a in starts:
        run_protocol(look_dry, LAD, a, eps, hs=HS, kinds=KINDS, **V121)
    todo = sorted(i for i in touched if i not in flag)
    print(f"first pass: {len(touched)} settings requested, {len(todo)} to flag", flush=True)
    with Pool(NPROC) as pool:
        for idx, val in pool.imap_unordered(_run, todo, chunksize=4):
            _record(idx, val)
    pickle.dump(flag, open(FLAGS, "wb"))
    print(f"flagged {len(flag)} runs in {time.time()-t0:.0f}s: "
          f"{sum(not v[2].get('converged', True) for v in flag.values())} did not converge, "
          f"{sum(not v[2].get('same', True) for v in flag.values())} differ from the stored cache",
          flush=True)

    RES = {}
    for mode, scr in [("released (default, no screening)", False)]:
        cert = miss = 0; pens, runs, calib, recs = [], [], [], set(); worst = 0.0
        for n, a in enumerate(starts):
            r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, **V121)
            runs.append(r["runs"])
            calib.append(r["protocol_cost"] / store[a][1])
            if r["certified"]:
                q = r["proposal"]; cert += 1; recs.add(q)
                e = abs(store[q][0][OB] - EX[OB]) / eps
                worst = max(worst, e); miss += e > 1
                pens.append(store[q][1] / c5)
            if n % 20 == 19:
                print(f"  {mode}: {n+1}/{len(starts)} starts done", flush=True)
        pickle.dump(flag, open(FLAGS, "wb"))
        seen = [i for i in store if abs(store[i][0][OB] - EX[OB]) <= eps]
        cbest = min(store[i][1] for i in seen)
        print(f"\n8-parameter L3, {mode} mode, energy 1e-2, {len(starts)} starts:")
        print(f"  certified {cert}, missed {miss}, worst error {worst:.2f} eps, "
              f"distinct recommendations {len(recs)}")
        print(f"  median runs {np.median(runs):.0f} (max {max(runs)}), median calibration cost "
              f"{np.median(calib):.0f} production runs at the start setting")
        print(f"  penalty against the 5-parameter grid optimum: median {np.median(pens):.2f}, "
              f"worst {max(pens):.2f}; against the cheapest attaining setting among the "
              f"{len(store)} 8-parameter runs made: median "
              f"{np.median([p * c5 / cbest for p in pens]):.2f}")
        RES[mode] = dict(starts=len(starts), accepted=cert, missed=int(miss), worst=worst,
                         distinct=len(recs), runs_med=float(np.median(runs)), runs_max=int(max(runs)),
                         calib_med=float(np.median(calib)),
                         pen5_med=float(np.median(pens)) if pens else None,
                         pen5_max=float(max(pens)) if pens else None,
                         pen8_med=float(np.median([p * c5 / cbest for p in pens])) if pens else None,
                         runs_seen=len(store),
                         unconverged_flagged=sum(not v[2].get("converged", True) for v in flag.values()))
    json.dump(RES, open(os.path.join(OUT, "eight_released.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
