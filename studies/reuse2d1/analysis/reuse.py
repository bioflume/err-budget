#!/usr/bin/env python3
"""The analysis of the reuse study (PROTOCOL.md 5-10), on the Mac:

    python3 analysis/reuse.py --runs results/<pack>/runs --out results/<pack>/res
                              [--starts-every K] [--engine DIR]

  1. accounting of the grid, the references and the check;
  2. the references, their uncertainty u and the targets graded (5, 6);
  3. the hardest member (7);
  4. (a): errbudget from the default start at the hardest member, its setting reused at
     every member, against the tutorial (b) and the cautious (c) setting (8, 9: R1);
  5. payback: S, K, mu, n* and the study's N* (9: R2), in CPU seconds and cells x iterations;
  6. (a, every start): from every converged start at the hardest member (--starts-every K:
     every K-th start only, for tests);
  7. (a'): the default-start search at every other member, reused at every member.
Writes reuse_summary.json.  Needs the errbudget engine: the article's src/ (default
~/OpenResearch/knobs-paper/knobs-paper/src, or $ERRBUDGET_SRC, or --engine).
"""
from __future__ import annotations

import argparse
import collections
import glob
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from cyl import study  # noqa: E402

COSTS = ("cpu", "cells_iters")


def engine(path=None):
    src = path or os.environ.get("ERRBUDGET_SRC") or os.path.expanduser("~/OpenResearch/knobs-paper/knobs-paper/src")
    sys.path.insert(0, src)
    import errbudget
    return errbudget


def load_runs(runs_dir):
    out = {}
    for p in glob.glob(os.path.join(runs_dir, "*.json")):
        name = os.path.basename(p)[:-5]
        if ".tmp" in name:
            continue
        with open(p) as f:
            out[name] = json.load(f)
    return out


def usable(x):
    return bool(x and x.get("status") == "ok" and x.get("converged") and math.isfinite(float(x.get("cd", math.nan))))


def fmt(x, nd=4):
    if x is None:
        return "-"
    if isinstance(x, float):
        if not math.isfinite(x):
            return str(x)
        return f"{x:.{nd}g}"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(fmt(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


# ------------------------------------------------------------------ 1. accounting
def accounting(runs):
    rows, total = [], {}
    for name, tasks in (("check", study.tasks_check()), ("grid", study.tasks_grid()), ("ref", study.tasks_ref())):
        rr = [r for t in tasks for r in t.runs]
        c = collections.Counter()
        for r in rr:
            x = runs.get(r.tag)
            c["stored"] += x is not None
            if x is None:
                c["missing"] += 1
            elif x["status"] == "failed":
                c["failed"] += 1
            elif x["converged"]:
                c["converged"] += 1
            else:
                c["at the cap"] += 1
        cpu = sum(float(runs[r.tag].get("cpu") or 0.0) for r in rr if r.tag in runs) / 3600
        rows.append([name, len(rr), c["stored"], c["converged"], c["at the cap"], c["failed"], c["missing"], round(cpu, 2)])
        total[name] = dict(expected=len(rr), simplefoam_core_h=cpu,
                           **{k: c[k] for k in ("stored", "converged", "at the cap", "failed", "missing")})
    return rows, total


# ------------------------------------------------------------------ 2. references
def references(runs):
    out = {}
    for re_ in study.MEMBERS:
        q, conv, it = {}, {}, {}
        for r in study.reference_runs(re_):
            x = runs.get(r.tag)
            conv[r.n_arc] = usable(x)
            q[r.n_arc] = float(x["cd"]) if usable(x) else math.nan
            it[r.n_arc] = x.get("iterations") if x else None
        q96, q128, q192 = q[96], q[128], q[192]
        d1, d2 = q128 - q96, q192 - q128
        p = math.log(abs(d1) / abs(d2)) / math.log(1.5) if d1 and d2 and math.isfinite(d1 + d2) else math.nan
        monotone = bool(d1 * d2 > 0) if math.isfinite(d1 + d2) else None
        pp = min(max(p, 1.0), 2.0) if (monotone and math.isfinite(p)) else 1.0
        u = abs(d2) / (1.5 ** pp - 1.0) / abs(q192) if math.isfinite(d2) and conv[192] else math.nan
        has_ref = bool(conv[192]) and math.isfinite(u)
        out[re_] = dict(q=q, converged=conv, iterations=it, p=p, p_used=pp, monotone=monotone, u=u,
                        q_ref=q192 if conv[192] else math.nan, has_ref=has_ref)
        if re_ in study.BENCHMARK and math.isfinite(q192):
            out[re_]["benchmark"] = study.BENCHMARK[re_]
            out[re_]["benchmark_dev"] = abs(q192 - study.BENCHMARK[re_]) / study.BENCHMARK[re_]
    return out


# ------------------------------------------------------------------ the grid of one member
def tensor(runs, re_):
    T = {}
    for idx in study.GRID:
        x = runs.get(study.run_at(re_, study.setting_of(idx)).tag)
        if x is None:
            T[idx] = (math.nan, {"cpu": math.inf, "cells_iters": math.inf}, dict(converged=False, missing=True))
            continue
        ok = usable(x)
        T[idx] = (float(x["cd"]) if ok else math.nan,
                  {"cpu": float(x.get("cpu") or 0.0), "cells_iters": float(x.get("cells_iters") or 0.0)},
                  dict(converged=ok, status=x["status"], iterations=x.get("iterations")))
    return T


def hardest_member(Ts):
    """PROTOCOL.md 7: the relative change of c_D from n_arc 24 to 32 at the default start."""
    s0 = study.START_IDX
    s1 = (s0[0] + 1,) + s0[1:]
    rows, worst = [], None
    for re_, T in Ts.items():
        a, b = T[s0], T[s1]
        conv = a[2]["converged"] and b[2]["converged"]
        ch = abs(b[0] - a[0]) / abs(a[0]) if conv else math.inf
        rows.append([re_, a[0], b[0], ch, conv])
        key = (ch, re_)
        if worst is None or key > worst[0]:
            worst = (key, re_)
    return worst[1], rows


def params_of(eb):
    hs = {"n_arc": lambda v: 1.0 / v, "n_coarsest": lambda v: 1.0 / v}
    return [eb.Parameter(n, v, kind=k, **({"h": hs[n]} if n in hs else {}))
            for n, v, k in zip(study.NAMES, study.LADDERS, study.KINDS)]


def search(eb, T, start_idx, target):
    """errbudget as a user calls it (tune, relative target) on the stored runs of one member."""
    used = []

    def solver(**kw):
        idx = study.index_of(kw)
        v, cost, info = T[idx]
        if not math.isfinite(cost["cpu"]):
            raise eb.RunFailed(f"no stored run for {kw}")
        if idx not in used:
            used.append(idx)
        return {"cd": v if math.isfinite(v) else 0.0, "converged": bool(info["converged"])}, cost["cpu"]

    res = eb.tune(solver, params_of(eb), target=target, observable="cd", start=study.setting_of(start_idx),
                  relative=True)
    out = dict(certified=bool(res.certified), error_estimate=float(res.error_estimate), eps_abs=float(res.target),
               limiting=res.limiting, runs=len(used), used=[list(i) for i in used],
               search={c: float(sum(T[i][1][c] for i in used if math.isfinite(T[i][1][c]))) for c in COSTS})
    if res.certified:
        out["setting"] = {k: (int(v) if k in ("n_arc", "n_coarsest") else float(v)) for k, v in res.settings.items()}
        out["index"] = list(study.index_of(res.settings))
    return out


def reuse(Ts, refs, idx, eps):
    """The setting idx at every member with a reference: (met, error / eps, converged)."""
    out = {}
    for re_, T in Ts.items():
        if not refs[re_]["has_ref"]:
            continue
        v, cost, info = T[tuple(idx)]
        qr = refs[re_]["q_ref"]
        err = abs(v - qr) / abs(qr) if math.isfinite(v) else math.inf
        out[re_] = dict(met=bool(info["converged"] and err <= eps), ratio=err / eps, cd=v,
                        converged=bool(info["converged"]), cost=cost)
    return out


def payback(Ts, h, idx_a, S):
    """PROTOCOL.md 9, R2, per cost measure."""
    idx_c = study.index_of(study.CAUTIOUS)
    out = {}
    for c in COSTS:
        ca_h, cc_h = Ts[h][tuple(idx_a)][1][c], Ts[h][idx_c][1][c]
        K = S[c] / ca_h
        mu = cc_h / ca_h
        save = [T[idx_c][1][c] - T[tuple(idx_a)][1][c] for T in Ts.values()]
        mean_save = float(np.mean(save))
        n_star = K / (mu - 1.0) if mu > 1 else math.inf
        N_star = S[c] / mean_save if mean_save > 0 else math.inf
        cum, paid_at = 0.0, None
        for k, s in enumerate(save, 1):
            cum += s
            if paid_at is None and cum >= S[c]:
                paid_at = k
        out[c] = dict(S=S[c], cost_a_hardest=ca_h, cost_c_hardest=cc_h, K=K, mu=mu, n_star=n_star,
                      mean_saving=mean_save, N_star=N_star, pays_within_study=bool(N_star <= len(Ts)),
                      paid_after_members_in_Re_order=paid_at,
                      mu_by_member=[T[idx_c][1][c] / T[tuple(idx_a)][1][c] for T in Ts.values()])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--out")
    ap.add_argument("--engine")
    ap.add_argument("--starts-every", type=int, default=1)
    a = ap.parse_args(argv)
    eb = engine(a.engine)
    runs = load_runs(a.runs)
    out_dir = a.out or os.path.join(os.path.dirname(os.path.abspath(a.runs)), "res")
    os.makedirs(out_dir, exist_ok=True)
    summary = dict(engine=getattr(eb, "__version__", "?"), engine_path=os.path.dirname(eb.__file__))
    print(f"errbudget {summary['engine']} from {summary['engine_path']}")

    print("\n== 1. Accounting")
    rows, summary["accounting"] = accounting(runs)
    print(md_table(["sweep", "runs", "stored", "converged", "at the cap", "failed", "missing",
                    "simpleFoam core-h"], rows))

    print("\n== 2. References (PROTOCOL.md 5)")
    refs = references(runs)
    summary["references"] = {str(k): v for k, v in refs.items()}
    print(md_table(["Re", "c_D n96", "n128", "n192", "iterations n192", "p", "p'", "u", "reference"],
                   [[re_, r["q"][96], r["q"][128], r["q"][192], r["iterations"][192], r["p"], r["p_used"],
                     r["u"], "yes" if r["has_ref"] else "NO"] for re_, r in refs.items()]))
    if 20.0 in refs and "benchmark_dev" in refs[20.0]:
        print(f"Re 20: n192 {refs[20.0]['q'][192]:.9f} against the benchmark {study.BENCHMARK[20.0]}: "
              f"{refs[20.0]['benchmark_dev']:.2e} relative")
    with_ref = [re_ for re_ in study.MEMBERS if refs[re_]["has_ref"]]
    u_max = max((refs[re_]["u"] for re_ in with_ref), default=math.nan)
    graded = {tg: bool(with_ref and tg >= 10 * u_max) for tg in study.TARGETS}
    summary["graded"] = {str(k): v for k, v in graded.items()}
    summary["u_max"] = u_max
    print(f"members with a reference: {len(with_ref)} of {len(study.MEMBERS)}; largest u {u_max:.3g}; "
          + "; ".join(f"target {tg:.1%}: {'GRADED' if g else 'NOT graded'}" for tg, g in graded.items()))

    Ts = {re_: tensor(runs, re_) for re_ in study.MEMBERS}
    print("\n== 3. The hardest member (PROTOCOL.md 7)")
    h, hrows = hardest_member(Ts)
    summary["hardest"] = h
    print(md_table(["Re", "c_D n24 (default start)", "c_D n32", "relative change", "both converged"], hrows))
    print(f"hardest member: Re = {h:g}")

    idx_b, idx_c = study.index_of(study.TUTORIAL), study.index_of(study.CAUTIOUS)
    summary["targets"] = {}
    for tg in study.TARGETS:
        print(f"\n== 4-5. Target {tg:.1%} ({'graded' if graded[tg] else 'NOT graded: information only'})")
        res = dict(graded=graded[tg])
        a_ = search(eb, Ts[h], study.START_IDX, tg)
        res["a"] = a_
        if not a_["certified"]:
            print(f"(a) DECLINED from the default start at Re {h:g} (estimate {a_['error_estimate']:.3g}, "
                  f"limiting {a_['limiting']}), {a_['runs']} runs")
        else:
            print(f"(a) at Re {h:g}: {a_['setting']} after {a_['runs']} runs "
                  f"(search {a_['search']['cpu']:.0f} CPU s); limiting {a_['limiting']}")
        rows = []
        settings = [("(a)", a_.get("index")), ("(b) tutorial", list(idx_b)), ("(c) cautious", list(idx_c))]
        reuses = {}
        for name, idx in settings:
            if idx is None:
                continue
            ru = reuse(Ts, refs, idx, tg)
            reuses[name] = ru
            met = sum(v["met"] for v in ru.values())
            worst = max((v["ratio"] for v in ru.values()), default=math.nan)
            rows.append([name, f"{study.setting_of(idx)}", f"{met}/{len(ru)}", worst,
                         ", ".join(f"{k:g}" for k, v in ru.items() if not v["met"]) or "-"])
        print(md_table(["setting", "values", "members met", "worst error / eps", "missed at Re"], rows))
        res["reuse"] = {k: {str(m): v for m, v in ru.items()} for k, ru in reuses.items()}
        if a_["certified"]:
            pb = payback(Ts, h, a_["index"], a_["search"])
            res["payback"] = pb
            for c in COSTS:
                p = pb[c]
                print(f"payback ({c}): S {p['S']:.4g}, K {p['K']:.3g}, mu {p['mu']:.3g} at Re {h:g} "
                      f"(median over members {np.median(p['mu_by_member']):.3g}), n* = K/(mu-1) {p['n_star']:.3g}; "
                      f"study: N* {p['N_star']:.3g} members ({'pays within the 13' if p['pays_within_study'] else 'does NOT pay within the 13'}; "
                      f"in Re order after {p['paid_after_members_in_Re_order']})")

        print(f"-- (a, every start) at Re {h:g}")
        starts = [i for i in study.GRID if Ts[h][i][2]["converged"]][::a.starts_every]
        every = collections.Counter()
        settings_seen, Ks, Ns, pairs_met, pairs = collections.Counter(), [], [], 0, 0
        for s in starts:
            r = search(eb, Ts[h], s, tg)
            if not r["certified"]:
                every["declined"] += 1
                continue
            every["accepted"] += 1
            settings_seen[tuple(r["index"])] += 1
            ru = reuse(Ts, refs, r["index"], tg)
            m = sum(v["met"] for v in ru.values())
            pairs_met += m
            pairs += len(ru)
            every["met at every member"] += m == len(ru)
            pb = payback(Ts, h, r["index"], r["search"])["cpu"]
            Ks.append(pb["K"])
            Ns.append(pb["N_star"])
        res["every_start"] = dict(starts=len(starts), **every, distinct_settings=len(settings_seen),
                                  pairs_met=pairs_met, pairs=pairs,
                                  K_median=float(np.median(Ks)) if Ks else math.nan,
                                  N_star_median=float(np.median(Ns)) if Ns else math.nan,
                                  N_star_share_within=float(np.mean([n <= len(Ts) for n in Ns])) if Ns else math.nan)
        e = res["every_start"]
        print(f"{e['starts']} converged starts: {e.get('accepted', 0)} accepted, {e.get('declined', 0)} declined, "
              f"{e['distinct_settings']} distinct settings; (start, member) pairs met {pairs_met}/{pairs}; "
              f"accepted searches whose setting meets the target at every member {e.get('met at every member', 0)}; "
              f"K median {e['K_median']:.3g}; N* median {e['N_star_median']:.3g} "
              f"({e['N_star_share_within']:.0%} within 13 members)")

        print("-- (a') the default-start search at each member, reused at every member (information)")
        rows, other = [], {}
        for re_ in study.MEMBERS:
            r = search(eb, Ts[re_], study.START_IDX, tg)
            if not r["certified"]:
                rows.append([re_, "declined", "-", "-"])
                other[str(re_)] = dict(certified=False)
                continue
            ru = reuse(Ts, refs, r["index"], tg)
            m = sum(v["met"] for v in ru.values())
            rows.append([re_, f"{r['setting']}", f"{m}/{len(ru)}",
                         ", ".join(f"{k:g}" for k, v in ru.items() if not v["met"]) or "-"])
            other[str(re_)] = dict(certified=True, setting=r["setting"], met=m, members=len(ru))
        print(md_table(["searched at Re", "setting", "members met", "missed at Re"], rows))
        res["other_members"] = other
        summary["targets"][str(tg)] = res

    path = os.path.join(out_dir, "reuse_summary.json")
    with open(path, "w") as f:
        json.dump(summary, f, indent=1, default=lambda o: None if isinstance(o, float) and not math.isfinite(o) else str(o))
    print(f"\nsummary written to {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
