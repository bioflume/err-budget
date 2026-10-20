"""Pooled results table: every case, observable and admitted target, one configuration.

    python pooled.py [config ...]        # configs defined in CONFIGS below

For each (case, observable, target) it replays the program from every usable
start and reports: starts; accepted; missed (accepted but outside the target);
declined although some usable, converged setting on the grid meets the target;
declined where none does; the median and worst cost penalty of the accepted
settings relative to the cheapest converged setting that meets the target; the
median number of runs; and how many accepted settings are runs whose nested
solver did not meet its tolerance.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from errbudget import run_protocol, V11, V12
try:
    from errbudget import V121
except ImportError:                      # the frozen program of the held-out tests
    V121 = dict(V12, joint_res=1.0)
import cases as C

# Program configurations compared in the article.  "v1.1" is the version of
# revision 2 (constant safety factor, unconverged runs accepted); "v1.2" is the
# default of this version.
CONFIGS = {
    # version 1.1 (revision 2): constant safety factor, unconverged runs accepted
    "v1.1": dict(fail_unconverged=False, kw=dict(V11)),
    "v1.1+gated": dict(fail_unconverged=False, kw=dict(V11, resid_margin=5.0)),
    # version 1.2, frozen for the held-out tests
    "v1.2": dict(fail_unconverged=True, kw=dict(V12)),
    # version 1.2.1, the default of the released program (amendment after H2)
    "v1.2.1": dict(fail_unconverged=True, kw=dict(V121)),
    # version 1.2 with one component removed (ablations of Table 3)
    "v1.2-accept-unconverged": dict(fail_unconverged=False, kw=dict(V12)),
    "v1.2-nolocal": dict(fail_unconverged=True, kw=dict(V12, local_order_check=False)),
    "v1.2-nojoint": dict(fail_unconverged=True, kw=dict(V12, joint_check=False)),
    "v1.2-credit": dict(fail_unconverged=True, kw=dict(V12, cancel_rule="credit")),
    "v1.2-constant": dict(fail_unconverged=True, kw=dict(V12, resid_margin=0.0)),
    "v1.2-nores": dict(fail_unconverged=True, kw=dict(V12, res_frac=0.0)),
    "v1.2-screened": dict(fail_unconverged=True, kw=dict(V12, use_screen=True)),
    # the nested-solver class removed: every nested-solver parameter declared as a
    # discretisation parameter (referee point 10)
    "v1.2-nested-as-disc": dict(fail_unconverged=True, kw=dict(V12), kinds_map={"nested-solver": "discretisation"}),
    # Revision 4: the same ablations taken relative to the released default, version
    # 1.2.1, so that Table 4 compares every row with the configuration of Table 3.
    "v1.2.1-accept-unconverged": dict(fail_unconverged=False, kw=dict(V121)),
    "v1.2.1-nolocal": dict(fail_unconverged=True, kw=dict(V121, local_order_check=False)),
    "v1.2.1-nojoint": dict(fail_unconverged=True, kw=dict(V121, joint_check=False)),
    "v1.2.1-credit": dict(fail_unconverged=True, kw=dict(V121, cancel_rule="credit")),
    "v1.2.1-constant": dict(fail_unconverged=True, kw=dict(V121, resid_margin=0.0)),
    "v1.2.1-nores": dict(fail_unconverged=True, kw=dict(V121, res_frac=0.0)),
    "v1.2.1-screened": dict(fail_unconverged=True, kw=dict(V121, use_screen=True)),
    "v1.2.1-nested-as-disc": dict(fail_unconverged=True, kw=dict(V121),
                                  kinds_map={"nested-solver": "discretisation"}),
    # Revision 4, referee point on GROMACS: version 1.2.1 without the round-off probe
    "v1.2.1-noprobe": dict(fail_unconverged=True, kw=dict(V121), no_floor=True),
    "v1.1-noprobe": dict(fail_unconverged=False, kw=dict(V11), no_floor=True),
}


def replay(case, ob, tg, cfg, starts=None):
    eps = tg * abs(case.exact[ob])
    fail_unc = cfg["fail_unconverged"]
    T = case.tensor

    def ok(i):
        return case.usable(i) and (case.converged(i) or not fail_unc)

    def look(i):
        i = tuple(i)
        v = T[i]
        if not np.isfinite(v[1]):
            return np.nan, np.inf
        if fail_unc and not case.converged(i):
            return np.nan, v[1]            # paid for, but treated as failed
        return v[0][ob], v[1]

    G = case.grid()
    err = {i: abs(T[i][0][ob] - case.exact[ob]) for i in G if case.usable(i)}
    # the optimum and reachability are judged over converged, usable runs, the
    # same for every configuration
    good = [i for i in G if case.usable(i) and case.converged(i)]
    feas = [i for i in good if err[i] <= eps]
    c_opt = min((T[i][1] for i in feas), default=np.nan)
    starts = [i for i in G if case.usable(i)] if starts is None else starts
    kw = dict(cfg["kw"])
    if case.floor is not None and not cfg.get("no_floor"):
        kw.setdefault("noise_floor", case.floor)
    kinds = cfg.get("kinds") or [cfg.get("kinds_map", {}).get(k, k) for k in case.KINDS]
    out = dict(case=case.key, ob=ob, target=tg, n=len(starts), reachable=bool(feas),
               accepted=0, missed=0, declined_reach=0, declined_unreach=0, start_failed=0,
               pen=[], runs=[], unconv=0, worst=0.0, misses=[], orders=[], hit_cap=0,
               max_steps=0, proposals={}, proposals_missed={}, unconv_from_converged=0,
               unconv_settings={}, reasons={}, feasible=len(feas), good=len(good))
    for a in starts:
        r = run_protocol(look, case.LAD, a, eps, hs=case.HS, kinds=kinds,
                         models=case.models, **kw)
        out["runs"].append(r["runs"])
        out["hit_cap"] += bool(r.get("hit_cap"))
        out["max_steps"] = max(out["max_steps"], int(r.get("iters", 0)))
        if r["certified"]:
            q = tuple(r["proposal"])
            out["accepted"] += 1
            # Revision 4: distinct accepted settings (counts over starts are not
            # independent samples; Section 5.3 promises the number of distinct settings)
            out["proposals"][q] = out["proposals"].get(q, 0) + 1
            e = err[q] / eps
            out["worst"] = max(out["worst"], e)
            if not case.converged(q):
                out["unconv"] += 1
                out["unconv_settings"][q] = out["unconv_settings"].get(q, 0) + 1
                out["unconv_from_converged"] += bool(case.converged(tuple(a)))
            out["orders"].append((_orders(r.get("budget"), case.KINDS), err[q] / eps))
            if e > 1:
                out["missed"] += 1
                out["misses"].append((a, q, e))
                out["proposals_missed"][q] = out["proposals_missed"].get(q, 0) + 1
            if feas:
                out["pen"].append(T[q][1] / c_opt)
        else:
            why = r.get("declined") or "not_reached"
            out["reasons"][why] = out["reasons"].get(why, 0) + 1
            if why == "start_failed":
                out["start_failed"] += 1
            elif feas:
                out["declined_reach"] += 1
            else:
                out["declined_unreach"] += 1
    return out


def _orders(lb, kinds):
    """Observed orders along the discretisation axes at the accepted setting: the
    order of the offset fit and that of the three-point local check (NaN where an
    axis was not fitted or the local changes alternate)."""
    if lb is None:
        return []
    out = []
    for k, kd in enumerate(kinds):
        if kd != "discretisation":
            continue
        pf = lb.pk[k] if k < len(getattr(lb, "pk", [])) else np.nan
        loc = getattr(lb, "local", {}).get(k)
        pl = loc["p"] if isinstance(loc, dict) else np.nan
        out.append((float(pf), float(pl)))
    return out


def summarise(o):
    pen = o["pen"]
    return dict(case=o["case"], ob=o["ob"], target=o["target"], n=o["n"],
                reachable=o["reachable"], accepted=o["accepted"], missed=o["missed"],
                declined_reach=o["declined_reach"], declined_unreach=o["declined_unreach"],
                start_failed=o["start_failed"],
                pen_med=float(np.median(pen)) if pen else float("nan"),
                pen_max=float(np.max(pen)) if pen else float("nan"),
                runs_med=float(np.median(o["runs"])) if o["runs"] else float("nan"),
                unconv=o["unconv"], worst=o["worst"], hit_cap=o["hit_cap"],
                max_steps=o["max_steps"],
                low_order=sum(1 for ords, e in o["orders"]
                              if any(np.isfinite(x) and x < 1 for pr in ords for x in pr)),
                low_order_missed=sum(1 for ords, e in o["orders"] if e > 1 and
                                     any(np.isfinite(x) and x < 1 for pr in ords for x in pr)),
                # revision 4 bookkeeping (no effect on any count above)
                distinct=len(o.get("proposals", {})),
                distinct_missed=len(o.get("proposals_missed", {})),
                accepted_settings=[[list(q), c] for q, c in sorted(o.get("proposals", {}).items())],
                unconv_distinct=len(o.get("unconv_settings", {})),
                unconv_from_converged=o.get("unconv_from_converged", 0),
                reasons=dict(o.get("reasons", {})),
                feasible=o.get("feasible"), good=o.get("good"))


def run_config(name, cases=None, verbose=True):
    cfg = CONFIGS[name]
    rows = []
    for case in (cases or C.all_cases()):
        for ob, tgs in case.targets.items():
            for tg in tgs:
                if not case.admitted(ob, tg):
                    continue
                t0 = time.time()
                o = replay(case, ob, tg, cfg)
                s = summarise(o)
                s["miss_list"] = [(list(a), list(q), e) for a, q, e in o["misses"][:20]]
                rows.append(s)
                if verbose:
                    print(f"  {case.key:<7} {ob:<16} {tg:<7.0e} n={s['n']:<4} acc={s['accepted']:<4} "
                          f"miss={s['missed']:<4} decl(reach)={s['declined_reach']:<4} "
                          f"startfail={s['start_failed']:<4} "
                          f"decl(unreach)={s['declined_unreach']:<4} pen={s['pen_med']:.2f}/"
                          f"{s['pen_max']:.2f} runs={s['runs_med']:.0f} unconv={s['unconv']} "
                          f"worst={s['worst']:.2f} ({time.time()-t0:.0f}s)", flush=True)
    return rows


if __name__ == "__main__":
    # python pooled.py CONFIG [CONFIG ...] [--only KEY,KEY] [--skip KEY,KEY]
    # With --only or --skip the rows of the selected cases replace those already in
    # data/pooled_<config>.json, so a case whose tensor is rebuilt can be redone alone.
    args = sys.argv[1:]
    only = skip = None
    if "--only" in args:
        i = args.index("--only"); only = set(args[i + 1].split(",")); del args[i:i + 2]
    if "--skip" in args:
        i = args.index("--skip"); skip = set(args[i + 1].split(",")); del args[i:i + 2]
    names = args or ["v1.2"]
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    # Revision 4: "2D-2h" is the refined DFG 2D-2 study of the Hamilton job J2; it is
    # loaded only when hamilton/results/j2_tensor.pkl exists, so the default pooled
    # table is unchanged until it does (select it with --only 2D-2h).
    loaders = [("L", C._manufactured), ("2D-1", lambda: [C._dfg1()]), ("2D-2", lambda: [C._dfg2()]),
               ("2D-2h", lambda: [C._dfg2_refined()]),
               ("LMP", C._lammps), ("GMX", C._gromacs), ("ABINIT", lambda: [C._abinit()]),
               ("QE", lambda: [C._qe()])]
    sel = lambda key: (only is None or key in only) and (skip is None or key not in skip)
    cases = []
    for _, load in loaders:
        try:
            got = load()
        except FileNotFoundError as exc:
            print(f"skipping: {exc}", flush=True)
            continue
        cases += [c for c in got if sel(c.key)]
    for name in names:
        print(f"=== configuration {name}: {CONFIGS[name]}", flush=True)
        rows = run_config(name, cases)
        path = os.path.join(out_dir, f"pooled_{name}.json")
        if (only or skip) and os.path.exists(path):
            done = {r["case"] for r in rows}
            old = [r for r in json.load(open(path)) if r["case"] not in done]
            rows = old + rows
        json.dump(rows, open(path, "w"), indent=1)
        tot = lambda k: sum(r[k] for r in rows)
        print(f"TOTAL {name}: replays {tot('n')}, accepted {tot('accepted')}, missed "
              f"{tot('missed')}, declined though reachable {tot('declined_reach')}, "
              f"declined unreachable {tot('declined_unreach')}, accepted unconverged "
              f"{tot('unconv')}", flush=True)
