"""Analysis of J2 (DFG 2D-2 refined, with relaxation time and window as parameters).

Reads $WORK/dfg2d2_hpc/*/series.npz, builds the 6-parameter tensor, writes
results/j2_tensor.pkl and prints:
  * the reference and its uncertainty (two other reference runs, and the
    relaxation/window variants of the finest one);
  * error against mesh level at the finest other settings, with observed orders;
  * the replay of the program (default configuration, and version 1.1) from
    every usable start at every admitted target.
    python analyse_j2.py [config ...]

Revision 4: the case specification is saved with the tensor (paper/cases.py reads
it as case "2D-2h"); the targets run from 1e-1 down (the admission rule decides
which are graded); error against the time step and the relaxation time and window;
the largest Courant numbers; the accounting of the settings (failed, not converged,
missing); configurations v1.2.1 (default), v1.2, v1.1, v1.2.1-credit and v1.2.1-nojoint.
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
# a smoke test must not overwrite the results of the real runs
OUT = os.path.join(HERE, "..", "results_smoke" if J.SMOKE else "results")
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


REF_STATUS = {}


def load():
    tensor, refs = {}, {}
    old = []
    LAD = [J.NA, J.DT, J.OUTER, J.PTOL, RELAX, WIN]
    for task in J.TASKS():
        d = os.path.join(RES, task["id"])
        f = os.path.join(d, "series.npz")
        res = json.load(open(os.path.join(d, "result.json"))) if os.path.exists(
            os.path.join(d, "result.json")) else None
        z = np.load(f) if os.path.exists(f) else None
        th = task["theta"]
        if res is not None and res.get("inlet") != "analytic":
            old.append(task["id"])
        if task["id"].startswith("ref_"):
            if res is not None:
                REF_STATUS[task["id"]] = dict(failed=res.get("failed"),
                                              nonconverged_steps=res.get("nonconverged_steps"),
                                              max_outer=res.get("max_outer"))
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
    if old and os.environ.get("ERRBUDGET_J2_ALLOW_OLD") != "1":
        raise SystemExit(f"{len(old)} results come from before the inlet fix of November 2026 "
                         f"(no inlet = analytic in result.json), e.g. {old[0]}: their mapped meshes "
                         "had the coarser mesh's inlet profile.  Move them aside and rerun "
                         "(ERRBUDGET_J2_ALLOW_OLD=1 analyses them anyway, for diagnosis only).")
    return LAD, tensor, refs


TARGETS = {"cd_max": [1e-1, 3e-2, 1e-2, 3e-3, 1e-3], "cl_max": [1e-1, 3e-2, 1e-2, 3e-3],
           "st": [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]}
CONFIGS = ["v1.2.1", "v1.2", "v1.1", "v1.2.1-credit", "v1.2.1-nojoint"]


def accounting(tensor):
    """Settings by outcome, and the base runs by status (for Table 2)."""
    missing, failed, diverged, unconv = [], 0, 0, 0
    for task in J.TASKS():
        if task["id"].startswith("ref_"):
            continue
        d = os.path.join(RES, task["id"])
        if not os.path.exists(os.path.join(d, "result.json")):
            missing.append(task["id"])
            continue
        r = json.load(open(os.path.join(d, "result.json")))
        failed += bool(r.get("failed"))
        diverged += bool(r.get("failed")) and "series_points" in r
        unconv += r.get("nonconverged_steps", 0) > 0
    usable = sum(np.isfinite(v[1]) for v in tensor.values())
    conv = sum(np.isfinite(v[1]) and v[2].get("converged") for v in tensor.values())
    print(f"base runs: {len(J.TASKS()) - len(J.REFS)} ({len(missing)} not finished), "
          f"{failed} stopped with an error ({diverged} of them after writing a partial series), "
          f"{unconv} with at least one step at the outer-iteration limit")
    print(f"settings (base run x relaxation x window): {len(tensor)}; produced a value {usable}; "
          f"converged {conv}; no value {len(tensor) - usable}")
    if missing:
        print("  NOT FINISHED (resubmit before using these results):", " ".join(missing[:20]),
              "..." if len(missing) > 20 else "")
    return dict(settings=len(tensor), usable=int(usable), converged=int(conv),
                no_value=int(len(tensor) - usable), base_runs_failed=failed,
                base_runs_diverged_after_series=diverged, base_runs_unconverged=unconv,
                missing=missing)


def courant_table():
    print("\nlargest Courant number of the runs (median and max over the tolerances), by mesh"
          " and time step; spin-ups:")
    for na in sorted(set(J.NA + [r["n_arc"] for r in J.REFS])):
        print(f"  spin-up n_arc {na}: {J.spun(na)}")
    out = {}
    for na in J.NA:
        row = []
        for dt in J.DT:
            co = []
            for o in J.OUTER:
                for p in J.PTOL:
                    f = os.path.join(RES, J.tid(dict(n_arc=na, dt=dt, outer_tol=o, p_tol=p)), "result.json")
                    if os.path.exists(f):
                        c = json.load(open(f)).get("courant_max")
                        if c is not None:
                            co.append(c)
            out[f"{na}|{dt:g}"] = co
            row.append(f"{np.median(co):6.2f}/{max(co):6.2f}" if co else "     -/-     ")
        print(f"  n_arc {na:>3}: " + "  ".join(f"dt {dt:g}: {r}" for dt, r in zip(J.DT, row)))
    return out


if __name__ == "__main__":
    LAD, tensor, refs = load()
    os.makedirs(OUT, exist_ok=True)
    ok = sum(np.isfinite(v[1]) for v in tensor.values())
    print(f"J2: {len(tensor)} settings ({ok} usable), {len(refs)} reference runs")
    acc = accounting(tensor)
    co = courant_table()
    if not refs:
        pickle.dump(dict(LAD=LAD, tensor=tensor, refs=refs), open(os.path.join(OUT, "j2_tensor.pkl"), "wb"))
        raise SystemExit("no reference runs yet")
    key0 = "ref_" + J.tid(J.REFS[0])
    RW = (RELAX[1], WIN[2])
    if key0 not in refs or refs[key0].get(RW) is None:
        pickle.dump(dict(LAD=LAD, tensor=tensor, refs=refs), open(os.path.join(OUT, "j2_tensor.pkl"), "wb"))
        raise SystemExit(f"the reference run {key0} is missing or has no complete shedding period "
                         f"in its window R {RW[0]:g}, W {RW[1]:g}")
    print("reference runs (failed / steps at the outer-iteration limit / largest outer count):")
    for k in sorted(REF_STATUS):
        st = REF_STATUS[k]
        print(f"  {k}: failed {bool(st['failed'])}, {st['nonconverged_steps']} steps at the limit, "
              f"max outer {st['max_outer']}")
    if REF_STATUS.get(key0, {}).get("nonconverged_steps"):
        raise SystemExit(f"the reference run {key0} has steps at the outer-iteration limit: "
                         "it is not converged; report it and do not grade against it")
    R0 = refs[key0][RW][0]
    var = []
    for k, v in refs.items():
        for rw, o in v.items():
            if o is not None and not (k == key0 and rw == RW):
                if rw == RW or k == key0:
                    var.append(o[0])
    unc = {o: max(abs(x[o] - R0[o]) / abs(R0[o]) for x in var) for o in R0}
    print("reference (n_arc 192, dt 1.25e-4, R 2, W 4):", {k: round(v, 6) for k, v in R0.items()})
    print("relative change of each other reference run and window from the reference:")
    for k, v in refs.items():
        for rw, o in sorted(v.items()):
            if o is not None and not (k == key0 and rw == RW) and (rw == RW or k == key0):
                print(f"  {k} R {rw[0]:g} W {rw[1]:g}: " + "  ".join(
                    f"{ob} {abs(o[0][ob] - R0[ob]) / abs(R0[ob]):.1e}" for ob in R0))
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
    print("\nrelative error against the time step at n_arc 96, tightest tolerances, R 2, W 4:")
    fin_dt = (len(LAD[0]) - 1,)
    for o in R0:
        e = [abs(tensor[fin_dt + (b, len(LAD[2]) - 1, len(LAD[3]) - 1, 1, 2)][0][o] - R0[o]) / abs(R0[o])
             for b in range(len(LAD[1]))]
        print(f"  {o:<7} " + "  ".join(f"{x:.2e}" for x in e))
        print(f"  {'order':<7} " + "  ".join(
            f"{np.log(e[j]/e[j+1])/np.log(LAD[1][j]/LAD[1][j+1]):8.2f}" for j in range(len(e) - 1)))
    print("\nrelative error against the relaxation time R (rows) and window W (columns) at the"
          " finest setting (n_arc 96, dt 2.5e-4, tightest tolerances):")
    for o in R0:
        print(f"  {o}")
        for r in range(len(LAD[4])):
            print(f"    R {LAD[4][r]:g}: " + "  ".join(
                f"W {LAD[5][w]:g} {abs(tensor[(len(LAD[0]) - 1, len(LAD[1]) - 1, len(LAD[2]) - 1, len(LAD[3]) - 1, r, w)][0][o] - R0[o]) / abs(R0[o]):.2e}"
                for w in range(len(LAD[5]))))
    HS = [[1.0 / n for n in LAD[0]], list(LAD[1]), list(LAD[2]), list(LAD[3]),
          list(LAD[4]), [1.0 / w for w in LAD[5]]]
    KINDS = ["discretisation"] * 2 + ["nested-solver"] * 2 + ["discretisation"] * 2
    MODELS = ["power", "power", "power", "power", "exponential", "power"]
    case = C.Case("2D-2h", "OpenFOAM DFG 2D-2 (refined)", tensor, LAD, HS, KINDS, dict(R0),
                  {k: list(v) for k, v in TARGETS.items()}, models=MODELS, unc=unc,
                  cost_note="core-seconds")
    spec = dict(LAD=LAD, HS=HS, KINDS=KINDS, exact=dict(R0), targets=TARGETS, models=MODELS,
                unc=unc, reference=key0, reference_window=RW)
    pickle.dump(dict(LAD=LAD, tensor=tensor, refs=refs, case=spec, R0=R0, unc=unc,
                     accounting=acc, courant=co),
                open(os.path.join(OUT, "j2_tensor.pkl"), "wb"))
    print("\nadmitted targets (at least ten times the reference uncertainty):",
          {o: [t for t in v if case.admitted(o, t)] for o, v in TARGETS.items()})
    # Every target is replayed; those below ten times the reference uncertainty are marked
    # (admitted False): they are graded against the reference only, as exp71_dfg2_detail.py
    # does for the 2D-2 study of revision 3, and the pooled table leaves them out.
    rows = []
    for cfg in (sys.argv[1:] or CONFIGS):
        for ob, tgs in case.targets.items():
            for tg in tgs:
                adm = case.admitted(ob, tg)
                s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
                s.update(config=cfg, admitted=adm)
                rows.append(s)
                print(f"  {cfg:<14} {ob:<7} {tg:<7.0e} starts {s['n']} start-failed {s['start_failed']} "
                      f"accepted {s['accepted']} ({s['distinct']} distinct) "
                      f"missed {s['missed']} declined(reachable) {s['declined_reach']} "
                      f"penalty {s['pen_med']:.2f}/{s['pen_max']:.2f} runs {s['runs_med']:.0f} "
                      f"worst {s['worst']:.2f}" + ("" if adm else "   [not admitted]"), flush=True)
    json.dump(dict(unc=unc, ref=R0, rows=rows, accounting=acc),
              open(os.path.join(OUT, "j2_replay.json"), "w"), indent=1)
