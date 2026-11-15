"""Replay of errbudget and of the procedures of practice on stored runs.

Copied from the article's repository (commit 4d32148) so that the airfoil study is
graded by the same code as the article's tables:

  seq, coord, ref_check   paper/baselines.py, verbatim (S, SJ and Ref. of Section 6.3)
  counting                paper/exp72_baselines_pooled.py, verbatim
  grade                   paper/exp72_baselines_pooled.py `cell`, with the decline
                          reasons and distinct settings of paper/pooled.py `replay`

Changes: `grade` takes the tensor, lists and reference value as arguments instead of a
`cases.Case`, accepts a subset of starts, and returns the per-start records as well as
the counts.  The rules are unchanged: a run whose solver did not converge is a failed
run for every procedure (NaN to the procedure, its cost still paid); a returned
setting reaches the target if its run converged and its error is within eps;
otherwise it misses; only the program can decline.
"""
from __future__ import annotations

import numpy as np


# ------------------------------------------------------------ paper/baselines.py
def seq(look, LAD, start, delta, order=None):
    """Sequential one-at-a-time convergence. Returns (proposal, runs)."""
    n = len(LAD)
    cur = list(start)
    runs, seen = 0, set()

    def ev(idx):
        nonlocal runs
        idx = tuple(idx)
        if idx not in seen:
            seen.add(idx)
            runs += 1
        return look(idx)[0]

    for k in (order or range(n)):
        q = ev(cur)
        while cur[k] + 1 < len(LAD[k]):
            trial = list(cur)
            trial[k] += 1
            q2 = ev(trial)
            if abs(q2 - q) <= delta:        # stopped moving: accept and freeze
                break
            cur, q = trial, q2
    return tuple(cur), runs


def coord(look, LAD, start, delta, max_pass=3):
    """SEQ plus a diagonal check, restarting if the diagonal step moves Q."""
    n = len(LAD)
    cur, runs = tuple(start), 0
    for _ in range(max_pass):
        cur, r = seq(look, LAD, cur, delta)
        runs += r
        diag = tuple(min(cur[k] + 1, len(LAD[k]) - 1) for k in range(n))
        if diag == cur:
            break
        runs += 2
        if abs(look(diag)[0] - look(cur)[0]) <= delta:
            break
        cur = diag
    return cur, runs


def ref_check(look, LAD, ALL, eps):
    """Converged-reference check: cheapest configuration within eps of the
    finest setting of every parameter. Costs the whole grid, which is what a
    reference-based study pays if it does not already own the sweep."""
    fine = tuple(len(l) - 1 for l in LAD)
    if not np.isfinite(look(fine)[0]):
        # the finest setting failed: a user would take the most refined one that ran
        fine = max(ALL, key=lambda i: (sum(i), look(i)[1]))
    qref = look(fine)[0]
    ok = [i for i in ALL if abs(look(i)[0] - qref) <= eps]
    if not ok:
        return fine, len(ALL)
    return min(ok, key=lambda i: look(i)[1]), len(ALL)


# ------------------------------------------------------------ paper/exp72_baselines_pooled.py
def counting(look):
    seen = {}

    def f(i):
        i = tuple(i)
        if i not in seen:
            seen[i] = look(i)
        return seen[i]
    f.seen = seen
    return f


PROCS = {"program": None, "seq-1": 1.0, "seq-3": 3.0, "seq-10": 10.0,
         "coord-1": 1.0, "coord-3": 3.0, "coord-10": 10.0, "ref": None}


def grade(tensor, LAD, HS, KINDS, qref, eps, run_protocol, cfg, starts=None, models=None, procs=PROCS):
    """Every procedure from every start (default: every setting whose run converged).

    tensor: idx -> (value, cost, {"converged": bool}); a NaN value or infinite cost is a
    run that cannot be used.  qref: the reference value; eps: the absolute target."""
    G = list(tensor)
    good = [i for i in G if np.isfinite(tensor[i][1]) and tensor[i][2].get("converged", True)
            and np.isfinite(tensor[i][0])]
    err = {i: abs(tensor[i][0] - qref) for i in good}
    feas = [i for i in good if err[i] <= eps]
    c_opt = min((tensor[i][1] for i in feas), default=np.nan)

    def look(i):
        i = tuple(i)
        v = tensor[i]
        if not (np.isfinite(v[1]) and v[2].get("converged", True) and np.isfinite(v[0])):
            return np.nan, (v[1] if np.isfinite(v[1]) else np.inf)
        return v[0], v[1]

    starts = good if starts is None else [s for s in starts if s in set(good)]
    out = {p: dict(reached=0, missed=0, declined=0, declined_reach=0, declined_unreach=0, pen=[], runs=[],
                   K=[], misses=[], proposals={}, reasons={}, records=[]) for p in procs}
    for a in starts:
        for p, div in procs.items():
            lk = counting(look)
            o = out[p]
            if p == "program":
                r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KINDS, models=models, **cfg)
                if not r["certified"]:
                    why = r.get("declined") or "not_reached"
                    o["declined"] += 1
                    o["declined_reach" if feas else "declined_unreach"] += 1
                    o["reasons"][why] = o["reasons"].get(why, 0) + 1
                    o["runs"].append(len(lk.seen))
                    o["records"].append(dict(start=a, returned=None, declined=why, runs=len(lk.seen)))
                    continue
                q = tuple(r["proposal"])
            elif p.startswith("seq"):
                q, _ = seq(lk, LAD, a, eps / div)
            elif p.startswith("coord"):
                q, _ = coord(lk, LAD, a, eps / div)
            else:
                q, _ = ref_check(lk, LAD, good, eps)
                for i in good:
                    lk(i)
            q = tuple(q)
            o["runs"].append(len(lk.seen))
            spent = sum(v[1] for v in lk.seen.values() if np.isfinite(v[1]))
            o["proposals"][q] = o["proposals"].get(q, 0) + 1
            rec = dict(start=a, returned=q, runs=len(lk.seen), spent=spent,
                       error=(err[q] / eps if q in err else None))
            if q in err and err[q] <= eps:
                o["reached"] += 1
                o["pen"].append(tensor[q][1] / c_opt)
                o["K"].append(spent / tensor[q][1])
                rec["K"] = spent / tensor[q][1]
            else:
                o["missed"] += 1
                o["misses"].append((a, q, (err[q] / eps) if q in err else None))
            o["records"].append(rec)
    return dict(procs=out, n_starts=len(starts), feasible=len(feas), good=len(good), c_opt=c_opt,
                optimum=(min(feas, key=lambda i: tensor[i][1]) if feas else None))


def summarise(res):
    rows = {}
    for p, o in res["procs"].items():
        rows[p] = dict(reached=o["reached"], missed=o["missed"], declined=o["declined"],
                       declined_reachable=o["declined_reach"], declined_unreachable=o["declined_unreach"],
                       penalty_median=float(np.median(o["pen"])) if o["pen"] else float("nan"),
                       penalty_worst=float(np.max(o["pen"])) if o["pen"] else float("nan"),
                       runs_median=float(np.median(o["runs"])) if o["runs"] else float("nan"),
                       K_median=float(np.median(o["K"])) if o["K"] else float("nan"),
                       distinct_settings=len(o["proposals"]), reasons=dict(o["reasons"]),
                       worst_miss=max((m[2] for m in o["misses"] if m[2] is not None), default=None))
    return rows
