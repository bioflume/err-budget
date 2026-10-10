"""Analysis of J2 (DFG 2D-2 refined, with relaxation time and window as parameters).

Reads $WORK/dfg2d2_hpc/*/series.npz, builds the 6-parameter tensor, writes
results/j2_tensor.pkl and prints:
  * the reference and its uncertainty (two other reference runs, and the
    relaxation/window variants of the finest one);
  * error against mesh level at the finest other settings, with observed orders;
  * the replay of the program (default configuration, and version 1.1) from
    every usable start at every admitted target.
    python analyse_j2.py [config ...]
"""
import itertools
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import dfg2d2_hpc as J
from openfoam_dfg2d2 import _peak, _frequency, FEATFLOW
import cases as C
import pooled as P

D, UBAR = 0.1, 1.0
RELAX, WIN = J.RELAX, J.WIN
OUT = os.path.join(HERE, "..", "results")
RES = os.path.join(os.environ.get("WORK", "."), "dfg2d2_hpc")


def observe(z, t0, R, W):
    t, cd, cl, ts, cpu, nc = z["t"], z["cd"], z["cl"], z["ts"], z["cpu"], z["nc"]
    a, b = t0 + R, t0 + R + W
    if len(t) == 0 or t[-1] < b - 1e-9:
        return None
    m = (t >= a - 1e-9) & (t <= b + 1e-9)
    tt, d, l = t[m], cd[m], cl[m]
    cdm, n1 = _peak(tt, d, l)
    clm, n2 = _peak(tt, l, l)
    st = _frequency(tt, l) * D / UBAR
    k = np.searchsorted(ts, b - 1e-9)
    k = min(k, len(ts) - 1)
    cost = float(cpu[k])
    conv = not bool(np.any(nc[: k + 1]))
    obs = dict(cd_max=cdm, cl_max=clm, st=st)
    if not all(np.isfinite(list(obs.values()))) or min(n1, n2) < 1:
        return None
    return obs, cost, conv


def load():
    tensor, refs = {}, {}
    LAD = [J.NA, J.DT, J.OUTER, J.PTOL, RELAX, WIN]
    for task in J.TASKS():
        d = os.path.join(RES, task["id"])
        f = os.path.join(d, "series.npz")
        res = json.load(open(os.path.join(d, "result.json"))) if os.path.exists(
            os.path.join(d, "result.json")) else None
        z = np.load(f) if os.path.exists(f) else None
        th = task["theta"]
        if task["id"].startswith("ref_"):
            if z is not None:
                refs[task["id"]] = {(R, W): observe(z, res["t0"], R, W)
                                    for R in RELAX for W in WIN if R + W <= J.T_TOT + 1e-9}
            continue
        base = tuple(LAD[k].index(th[n]) for k, n in enumerate(
            ["n_arc", "dt", "outer_tol", "p_tol"]))
        for r, R in enumerate(RELAX):
            for w, W in enumerate(WIN):
                idx = base + (r, w)
                o = observe(z, res["t0"], R, W) if (z is not None and res) else None
                if o is None:
                    tensor[idx] = ({k: np.nan for k in ("cd_max", "cl_max", "st")}, np.inf,
                                   dict(converged=False))
                else:
                    tensor[idx] = (o[0], o[1], dict(converged=o[2]))
    return LAD, tensor, refs


if __name__ == "__main__":
    LAD, tensor, refs = load()
    os.makedirs(OUT, exist_ok=True)
    pickle.dump(dict(LAD=LAD, tensor=tensor, refs=refs), open(os.path.join(OUT, "j2_tensor.pkl"), "wb"))
    ok = sum(np.isfinite(v[1]) for v in tensor.values())
    print(f"J2: {len(tensor)} settings ({ok} usable), {len(refs)} reference runs")
    if not refs:
        raise SystemExit("no reference runs yet")
    key0 = "ref_" + J.tid(J.REFS[0])
    RW = (RELAX[1], WIN[2])
    R0 = refs[key0][RW][0]
    var = []
    for k, v in refs.items():
        for rw, o in v.items():
            if o is not None and not (k == key0 and rw == RW):
                if rw == RW or k == key0:
                    var.append(o[0])
    unc = {o: max(abs(x[o] - R0[o]) / abs(R0[o]) for x in var) for o in R0}
    print("reference (n_arc 192, dt 1.25e-4, R 2, W 4):", {k: round(v, 6) for k, v in R0.items()})
    print("relative uncertainty (other references and windows):", {k: f"{v:.1e}" for k, v in unc.items()})
    print("FEATFLOW level 6:", FEATFLOW, " differences:",
          {o: f"{abs(R0[o]-FEATFLOW[o])/FEATFLOW[o]:.1e}" for o in R0})
    fin = (len(LAD[1]) - 1, len(LAD[2]) - 1, len(LAD[3]) - 1, 1, 2)
    print("\nrelative error against n_arc at dt 2.5e-4, tightest tolerances, R 2, W 4:")
    for o in R0:
        e = [abs(tensor[(a,) + fin][0][o] - R0[o]) / abs(R0[o]) for a in range(len(LAD[0]))]
        print(f"  {o:<7} " + "  ".join(f"{x:.2e}" for x in e))
        print(f"  {'order':<7} " + "  ".join(
            f"{np.log(e[j]/e[j+1])/np.log(LAD[0][j+1]/LAD[0][j]):8.2f}" for j in range(len(e) - 1)))
    case = C.Case("2D-2h", "OpenFOAM DFG 2D-2 (refined)", tensor, LAD,
                  [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), list(LAD[3]),
                   list(LAD[4]), [1.0 / w for w in LAD[5]]],
                  ["discretisation"] * 2 + ["nested-solver"] * 2 + ["discretisation"] * 2,
                  dict(R0), {"cd_max": [1e-2, 3e-3, 1e-3], "cl_max": [3e-2, 1e-2, 3e-3],
                             "st": [1e-2, 3e-3, 1e-3]},
                  models=["power", "power", "power", "power", "exponential", "power"],
                  unc=unc, cost_note="core-seconds")
    rows = []
    for cfg in (sys.argv[1:] or ["v1.2", "v1.1"]):
        for ob, tgs in case.targets.items():
            for tg in tgs:
                if not case.admitted(ob, tg):
                    print(f"  {cfg} {ob} {tg:.0e}: below ten times the reference uncertainty")
                    continue
                s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
                s["config"] = cfg
                rows.append(s)
                print(f"  {cfg:<6} {ob:<7} {tg:<7.0e} starts {s['n']} accepted {s['accepted']} "
                      f"missed {s['missed']} declined(reachable) {s['declined_reach']} "
                      f"penalty {s['pen_med']:.2f}/{s['pen_max']:.2f} runs {s['runs_med']:.0f} "
                      f"worst {s['worst']:.2f}", flush=True)
    json.dump(dict(unc=unc, ref=R0, rows=rows), open(os.path.join(OUT, "j2_replay.json"), "w"),
              indent=1)
