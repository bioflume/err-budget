"""Analysis of J1 (DFG 2D-1 with a current OpenFOAM release).

    python analyse_j1.py [config ...]          # default: v1.2.1 v1.1

Prints, and writes to ../results/:
  j1_fo.json        forceCoeffs function object against the article's post-processing
                    of the fields, n_arc 8-192 at the tightest settings, with the wall
                    distance of the two cell centres and the viscous drag from each;
  j1_stop.json      why the linear solves stopped in the 840 grid runs (relTol 0.1),
                    by the ratio of the absolute tolerance to residualControl;
  j1_reltol.json    the relTol study at n_arc 32: relTol 0.1 (grid), 0.01 and 0;
  dfg_tensor_current.pkl, dfg_tensor_current_fo.pkl
                    the grid with the current release (field post-processing, and
                    c_D, c_L from the function object);
  j1_replay.json    the program on the grid with the v1912 tensor of the article and
                    with the two current-release tensors.
"""
import json
import os
import pickle
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import dfg2d1_hpc as J
from openfoam_cylinder import REF
import cases as C
import pooled as P

RES = os.path.join(os.environ.get("WORK", "."), "dfg2d1_hpc")
# a smoke test must not overwrite the results of the real runs
OUT = os.path.join(HERE, "..", "results_smoke" if J.SMOKE else "results")
os.makedirs(OUT, exist_ok=True)
CLASSES = ["abs", "rel", "both", "zero", "max", "other"]


def _load(tid):
    f = os.path.join(RES, tid, "result.json")
    return json.load(open(f)) if os.path.exists(f) else None


def fo_table():
    print("forces function object against the article's post-processing (current release)")
    print(f"{'n_arc':>6} {'iters':>6} {'conv':>5} {'Cd (FO)':>12} {'Cd (fields)':>12} {'Cd (v1.1)':>12}"
          f" {'err FO':>9} {'err fields':>10} {'err v1.1':>9} {'Cl FO':>10} {'Cl fields':>10}")
    fo = {}
    for n in J.FO_LEVELS:
        r = _load(f"fo_na{n}")
        if r is None:
            g = os.path.join(RES, f"fo_na{n}", "wall_geometry.json")
            if os.path.exists(g):
                fo[n] = dict(wall=json.load(open(g)), missing="solver run")
            continue
        fo[n] = r
        print(f"{n:>6} {r['iters']:>6} {str(r['converged'])[0]:>5} {r.get('cd_fo', np.nan):12.6f} "
              f"{r['cd_fields']:12.6f} {r['cd_fields_vertexmean']:12.6f} "
              f"{r.get('err_cd_fo', np.nan):9.2e} {r['err_cd_fields']:10.2e} "
              f"{r['err_cd_fields_vertexmean']:9.2e} {r.get('cl_fo', np.nan):10.6f} {r['cl_fields']:10.6f}")
    print("\nwall-cell centres: wall distance of the vertex mean relative to the centroid,"
          " and the viscous drag from each")
    print(f"{'n_arc':>6} {'d centroid':>11} {'offset min':>11} {'median':>9} {'max':>9}"
          f" {'cD visc (c)':>12} {'cD visc (vm)':>13} {'diff/cD':>9} {'n*offset':>9}")
    for n, r in sorted(fo.items()):
        w = r.get("wall", {})
        if not w:
            continue
        dv = (w.get("cd_viscous_vertexmean", np.nan) - w.get("cd_viscous_centroid", np.nan)) / REF["cd"]
        print(f"{n:>6} {w['d_centroid_median']:11.3e} {w['rel_offset_min']:11.3e} "
              f"{w['rel_offset_median']:9.3e} {w['rel_offset_max']:9.3e} "
              f"{w.get('cd_viscous_centroid', np.nan):12.6f} {w.get('cd_viscous_vertexmean', np.nan):13.6f} "
              f"{dv:9.2e} {n * w['rel_offset_median']:9.3f}")
    print("  (n*offset constant means the offset is a first-order fraction of the wall distance)")
    json.dump(fo, open(os.path.join(OUT, "j1_fo.json"), "w"), indent=1)
    return fo


def grid_tensors():
    cur, fo, raw = {}, {}, {}
    for t in J.TASKS():
        if t["kind"] != "grid":
            continue
        r = _load(t["id"])
        if r is None:
            continue
        i = tuple(r["idx"])
        info = dict(r["info"])
        c = r["cost"]
        if np.isfinite(c) and info.get("cells") and info.get("iters"):
            c = float(info["cells"] * info["iters"])      # as in the article (cases._dfg1)
        raw[i] = r
        obs = {k: r["obs"][k] for k in ("cd", "cl", "dp", "dp_cell", "cd_vm", "cl_vm")}
        cur[i] = (obs, c, info)
        ofo = dict(obs, cd=r["obs"].get("cd_fo", np.nan), cl=r["obs"].get("cl_fo", np.nan))
        fo[i] = (ofo, c if np.isfinite(ofo["cd"]) else np.inf, info)
    return cur, fo, raw


def stop_summary(raw):
    """Stopping reasons of the linear solves in the last 20 outer iterations of each grid
    run, pooled by log10(tolerance / residualControl) for p and for U."""
    agg = {"p": defaultdict(lambda: defaultdict(int)), "U": defaultdict(lambda: defaultdict(int))}
    runs = defaultdict(int)
    for i, r in raw.items():
        st = r["info"].get("stop")
        if not st:
            continue
        rc, t = r["theta"]["res_ctrl"], r["theta"]["lin_tol"]
        # p: tolerance t against residualControl rc; U: 10 t against rc / 10
        keys = {"p": int(round(np.log10(t / rc))), "U": int(round(np.log10(100 * t / rc)))}
        for g in ("p", "U"):
            runs[(g, keys[g])] += 1
            for c, n in st["last"][g].items():
                agg[g][keys[g]][c] += n
    out = {}
    print("\nwhy the linear solves stopped, last 20 outer iterations of the 840 grid runs"
          " (relTol 0.1); fractions of the solves")
    for g, lab in (("p", "GAMG (p): log10(tolerance / residualControl_p)"),
                   ("U", "smoothSolver (U): log10(tolerance_U / residualControl_U) = log10(100 t / rc)")):
        print(f"  {lab}")
        print(f"  {'ratio':>6} {'runs':>5} " + " ".join(f"{c:>7}" for c in CLASSES))
        out[g] = {}
        for k in sorted(agg[g]):
            tot = sum(agg[g][k].values())
            fr = {c: agg[g][k].get(c, 0) / tot for c in CLASSES} if tot else {}
            out[g][k] = dict(runs=runs[(g, k)], solves=tot, fractions=fr)
            print(f"  {k:>6} {runs[(g, k)]:>5} " + " ".join(f"{fr.get(c, 0):7.3f}" for c in CLASSES))
    print("  abs: stopped on the absolute tolerance only; rel: on relTol only; zero: no iteration"
          " (initial residual already below the tolerance); max: maxIter")
    json.dump(out, open(os.path.join(OUT, "j1_stop.json"), "w"), indent=1)
    return out


def reltol_study(raw):
    """relTol 0.1 (grid, n_arc 32, n_coarsest 10), 0.01 and 0 at every residualControl x
    GAMG tolerance; the change of the answer when relTol is removed."""
    ia = J.LAD[0].index(J.RT_NARC) if J.RT_NARC in J.LAD[0] else None
    ic = J.LAD[3].index(J.RT_NC)
    rows = []
    print(f"\nrelTol study at n_arc {J.RT_NARC}: relative errors and outer iterations")
    print(f"{'resCtrl':>8} {'GAMG tol':>9} {'relTol':>7} {'iters':>6} {'conv':>5} {'err cd':>9} "
          f"{'err cl':>9} {'err dp':>9} {'|d cd| vs 0.1':>14} {'|d cl|':>9} {'|d dp|':>9}")
    for rc in J.RT_RC:
        for t in J.RT_TOL:
            base = None
            if ia is not None:
                base = raw.get((ia, J.LAD[1].index(rc), J.LAD[2].index(t), ic))
            variants = [(0.1, base)] + [(rel, _load(f"rt_na{J.RT_NARC}_rc{rc:g}_t{t:g}_rel{rel:g}"))
                                        for rel in J.RT_REL]
            for rel, r in variants:
                if r is None:
                    continue
                o, info = r["obs"], r["info"]
                e = {k: abs(o[k] - REF[k]) / abs(REF[k]) for k in ("cd", "cl", "dp")}
                d = {k: (abs(o[k] - base["obs"][k]) / abs(REF[k]) if base else np.nan)
                     for k in ("cd", "cl", "dp")}
                row = dict(res_ctrl=rc, lin_tol=t, relTol=rel, iters=info.get("iters"),
                           converged=info.get("converged"), err=e, change_from_reltol_01=d,
                           stop_last=info.get("stop", {}).get("last"))
                rows.append(row)
                print(f"{rc:8.0e} {t:9.0e} {rel:7.2g} {str(info.get('iters')):>6} "
                      f"{str(info.get('converged'))[0]:>5} {e['cd']:9.2e} {e['cl']:9.2e} {e['dp']:9.2e} "
                      f"{d['cd']:14.2e} {d['cl']:9.2e} {d['dp']:9.2e}")
    json.dump(rows, open(os.path.join(OUT, "j1_reltol.json"), "w"), indent=1)
    return rows


def compare_post(cur, fo):
    both = [i for i in cur if np.isfinite(cur[i][1]) and np.isfinite(fo[i][0]["cd"])
            and cur[i][2].get("converged")]
    if not both:
        return
    dcd = np.array([abs(fo[i][0]["cd"] - cur[i][0]["cd"]) / REF["cd"] for i in both])
    dcl = np.array([abs(fo[i][0]["cl"] - cur[i][0]["cl"]) / REF["cl"] for i in both])
    print(f"\nfunction object against field post-processing over {len(both)} converged grid runs:"
          f" |d cD|/cD median {np.median(dcd):.2e} max {dcd.max():.2e};"
          f" |d cL|/cL median {np.median(dcl):.2e} max {dcl.max():.2e}")


def replay(cur, fo, configs):
    base = C._dfg1()
    mk = lambda key, lab, t: C.Case(key, lab, t, base.LAD, base.HS, base.KINDS, base.exact,
                                    base.targets, cost_note=base.cost_note)
    variants = [("v1912 (article)", base), ("current, fields", mk("2D-1c", "current release", cur))]
    missing = [i for i in cur if np.isfinite(cur[i][1]) and not np.isfinite(fo[i][0]["cd"])]
    if not missing:
        variants.append(("current, function object", mk("2D-1fo", "current release, FO", fo)))
    else:
        print(f"  function-object forces missing in {len(missing)} runs with fields "
              f"(e.g. {missing[:3]}): replay with function-object forces skipped")
    rows = []
    for cfg in configs:
        for lab, case in variants:
            for ob, tgs in case.targets.items():
                if lab.endswith("function object") and ob == "dp":
                    continue                       # the function object gives no pressure difference
                for tg in tgs:
                    s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
                    s.update(config=cfg, release=lab)
                    rows.append(s)
                    print(f"  {cfg:<7} {lab:<26} {ob:<3} {tg:.0e}: starts {s['n']} accepted "
                          f"{s['accepted']} ({s['distinct']} distinct) missed {s['missed']} "
                          f"declined {s['declined_reach']}+{s['declined_unreach']} "
                          f"penalty {s['pen_med']:.2f} runs {s['runs_med']:.0f}", flush=True)
    json.dump(rows, open(os.path.join(OUT, "j1_replay.json"), "w"), indent=1)


if __name__ == "__main__":
    fo_table()
    cur, fo, raw = grid_tensors()
    n = len(J.GRID)
    print(f"\ngrid with the current release: {len(cur)}/{n} settings, "
          f"{sum(1 for v in cur.values() if not np.isfinite(v[1]))} failed, "
          f"{sum(1 for v in cur.values() if np.isfinite(v[1]) and not v[2].get('converged'))} not converged")
    if raw:
        stop_summary(raw)
        reltol_study(raw)
    if len(cur) == n and not J.SMOKE:
        pickle.dump(cur, open(os.path.join(OUT, "dfg_tensor_current.pkl"), "wb"))
        pickle.dump(fo, open(os.path.join(OUT, "dfg_tensor_current_fo.pkl"), "wb"))
        compare_post(cur, fo)
        replay(cur, fo, sys.argv[1:] or ["v1.2.1", "v1.1"])
