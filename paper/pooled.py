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
from knobs import run_protocol, V11, V12
import cases as C

# Program configurations compared in the article.  "v1.1" is the version of
# revision 2 (constant safety factor, unconverged runs accepted); "v1.2" is the
# default of this version.
CONFIGS = {
    # version 1.1 (revision 2): constant safety factor, unconverged runs accepted
    "v1.1": dict(fail_unconverged=False, kw=dict(V11)),
    "v1.1+gated": dict(fail_unconverged=False, kw=dict(V11, resid_margin=5.0)),
    # version 1.2, the default of this article
    "v1.2": dict(fail_unconverged=True, kw=dict(V12)),
    # version 1.2 with one component removed (ablations of Table 3)
    "v1.2-accept-unconverged": dict(fail_unconverged=False, kw=dict(V12)),
    "v1.2-nolocal": dict(fail_unconverged=True, kw=dict(V12, local_order_check=False)),
    "v1.2-nojoint": dict(fail_unconverged=True, kw=dict(V12, joint_check=False)),
    "v1.2-credit": dict(fail_unconverged=True, kw=dict(V12, cancel_rule="credit")),
    "v1.2-constant": dict(fail_unconverged=True, kw=dict(V12, resid_margin=0.0)),
    "v1.2-nores": dict(fail_unconverged=True, kw=dict(V12, res_frac=0.0)),
    "v1.2-screened": dict(fail_unconverged=True, kw=dict(V12, use_screen=True)),
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
    if case.floor is not None:
        kw.setdefault("noise_floor", case.floor)
    out = dict(case=case.key, ob=ob, target=tg, n=len(starts), reachable=bool(feas),
               accepted=0, missed=0, declined_reach=0, declined_unreach=0, start_failed=0,
               pen=[], runs=[], unconv=0, worst=0.0, misses=[])
    for a in starts:
        r = run_protocol(look, case.LAD, a, eps, hs=case.HS, kinds=case.KINDS,
                         models=case.models, **kw)
        out["runs"].append(r["runs"])
        if r["certified"]:
            q = tuple(r["proposal"])
            out["accepted"] += 1
            e = err[q] / eps
            out["worst"] = max(out["worst"], e)
            out["unconv"] += not case.converged(q)
            if e > 1:
                out["missed"] += 1
                out["misses"].append((a, q, e))
            if feas:
                out["pen"].append(T[q][1] / c_opt)
        elif r.get("declined") == "start_failed":
            out["start_failed"] += 1
        elif feas:
            out["declined_reach"] += 1
        else:
            out["declined_unreach"] += 1
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
                unconv=o["unconv"], worst=o["worst"])


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
    names = sys.argv[1:] or ["v1.2"]
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    for name in names:
        print(f"=== configuration {name}: {CONFIGS[name]}", flush=True)
        rows = run_config(name)
        json.dump(rows, open(os.path.join(out_dir, f"pooled_{name}.json"), "w"), indent=1)
        tot = lambda k: sum(r[k] for r in rows)
        print(f"TOTAL {name}: replays {tot('n')}, accepted {tot('accepted')}, missed "
              f"{tot('missed')}, declined though reachable {tot('declined_reach')}, "
              f"declined unreachable {tot('declined_unreach')}, accepted unconverged "
              f"{tot('unconv')}", flush=True)
