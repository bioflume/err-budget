"""Replay of the program and of the procedures of baselines.py from every starting setting
of a stored study; used by the builders exp48_cavity.py, exp53_qe.py and exp56_dfg2d2.py.
"""
import itertools

import numpy as np

from errbudget import run_protocol
from baselines import seq, coord, ref_check


def grid(LAD):
    return list(itertools.product(*[range(len(l)) for l in LAD]))


def method_cell(tensor, LAD, HS, KINDS, ob, QS, eps, starts=None, **kw):
    """Replay the method from every start; grade against the exact value QS."""
    ALL = [i for i in grid(LAD) if np.isfinite(tensor[i][1])]
    starts = ALL if starts is None else starts
    err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
    cost = {i: tensor[i][1] for i in ALL}
    feas = [i for i in ALL if err[i] <= eps]
    c_opt = min((cost[i] for i in feas), default=np.nan)
    look = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    out = dict(n=len(starts), attainable=bool(feas), certified=0, missed=0,
               declined=0, false_decl=0, pen=[], runs=[], calib=[], worst_err=0.0,
               recs=set(), misses=[])
    for a in starts:
        r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, **kw)
        out["runs"].append(r["runs"])
        out["calib"].append(r["protocol_cost"] / cost[a])
        if r["certified"]:
            q = r["proposal"]
            out["certified"] += 1
            out["recs"].add(q)
            out["worst_err"] = max(out["worst_err"], err[q] / eps)
            if err[q] > eps:
                out["missed"] += 1
                out["misses"].append((a, q, err[q] / eps, abs(r["budget"].E0) / eps))
            if feas:
                out["pen"].append(cost[q] / c_opt)
        else:
            out["declined"] += 1
            out["false_decl"] += bool(feas)
    out["c_opt"] = c_opt
    out["opt"] = min(feas, key=lambda i: cost[i]) if feas else None
    return out


def baseline_cell(tensor, LAD, ob, QS, eps, starts=None, deltas=(1.0, 0.1)):
    """Sequential convergence and coordinate search at delta = eps and eps/10,
    and the converged-reference check, from the same starts."""
    ALL = [i for i in grid(LAD) if np.isfinite(tensor[i][1])]
    starts = ALL if starts is None else starts
    look = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1]) \
        if np.isfinite(tensor[tuple(i)][1]) else (np.nan, np.inf)
    err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
    feas = [i for i in ALL if err[i] <= eps]
    c_opt = min((tensor[i][1] for i in feas), default=np.nan)
    res = {}
    for name, fn in [("sequential", seq), ("coordinate search", coord)]:
        for dl in deltas:
            ok, pen, runs = 0, [], []
            for a in starts:
                q, n = fn(look, LAD, a, dl * eps)
                runs.append(n)
                if not np.isfinite(tensor[q][1]):
                    continue
                ok += err[q] <= eps
                pen.append(tensor[q][1] / c_opt)
            lab = f"{name} @eps" + ("" if dl == 1.0 else f"/{1/dl:.0f}")
            res[lab] = dict(attained=ok, n=len(starts), pen=float(np.median(pen)),
                            runs=float(np.median(runs)))
    q, n = ref_check(look, LAD, ALL, eps)
    res["reference acceptance"] = dict(attained=int(err[q] <= eps), n=1,
                                       pen=tensor[q][1] / c_opt, runs=n)
    return res


def fmt_method(o):
    pen = np.array(o["pen"]) if o["pen"] else np.array([np.nan])
    return (f"certified {o['certified']}/{o['n']}, missed {o['missed']}, "
            f"declined {o['declined']} ({o['false_decl']} attainable), "
            f"distinct recommendations {len(o['recs'])}, penalty median "
            f"{np.median(pen):.2f} worst {np.max(pen):.2f}, runs median "
            f"{np.median(o['runs']):.0f}, worst error {o['worst_err']:.2f} eps")
