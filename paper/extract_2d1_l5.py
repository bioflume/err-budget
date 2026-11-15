"""Extract the plotted values of fig_teaser, fig_sup_dfg, fig_sup_search and fig_l5_stall.

Data sources (all in the article's repository, PAPER_REPO, commit 4d32148, read only):
  paper/data/dfg_tensor_v12.pkl  OpenFOAM 2D-1 tensor, version-1.2 post-processing, loaded
                                 through paper/cases.py:_dfg1. "cd", "cl", "dp": forces from
                                 the exact wall-cell centroids and the second-order wall
                                 pressure; "cd_vm": drag with the vertex-mean cell centres.
  paper/data/dfg_fine.json       2D-1 at n_arc = 96, 128 (tightest settings, beyond the study).
  paper/data/cavity_tensor.pkl   L5, through paper/grids.py.
  src/errbudget                  the engine; the search is rerun with V121 exactly as
                                 paper/exp63_figures.py:fig_search does (v11 Fig. 4).
Writes figures/data_2d1_l5.json. Run with PYTHONDONTWRITEBYTECODE=1 so that nothing is
written into PAPER_REPO:
  PYTHONDONTWRITEBYTECODE=1 python3 figures/extract_2d1_l5.py
"""
import json
import os

import numpy as np

import cases as C                                   # noqa: E402
import grids                                        # noqa: E402
from errbudget import LocalBudget, run_protocol, V121   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def dfg_mesh():
    c = C._dfg1(cost="deterministic")
    T, LAD = c.tensor, c.LAD
    na = LAD[0]
    tight = lambda a: (a, len(LAD[1]) - 1, len(LAD[2]) - 1, 0)
    tut = lambda a: (a, 0, LAD[2].index(1e-6), 0)
    out = {"n_arc": na, "exact": c.exact, "tightest_setting": [LAD[1][-1], LAD[2][-1], LAD[3][0]],
           "tutorial_setting": [LAD[1][0], 1e-6, LAD[3][0]]}
    for ob in ["cd", "dp", "cl", "cd_vm"]:
        ex = c.exact[ob.split("_")[0]]
        out["tight_" + ob] = [abs(T[tight(a)][0][ob] - ex) / abs(ex) for a in range(len(na))]
        if ob != "cd_vm":
            out["tut_" + ob] = [abs(T[tut(a)][0][ob] - ex) / abs(ex) for a in range(len(na))]
    fj = json.load(open(os.path.join(C.D, "dfg_fine.json")))
    nf = [int(n) for n in sorted(fj, key=int) if int(n) > max(na)]
    out["n_arc_fine"] = nf
    for ob in ["cd", "dp", "cl"]:
        ex = c.exact[ob]
        out["fine_" + ob] = [abs(fj[str(n)]["obs"][ob] - ex) / abs(ex) for n in nf]
    return out


def search():
    """As exp63_figures.py:fig_search (v11 Fig. 4), unchanged."""
    c = C._dfg1(cost="deterministic")
    T, LAD, ob = c.tensor, c.LAD, "cl"
    eps = 3e-2 * abs(c.exact[ob])
    look = lambda i: ((T[tuple(i)][0][ob], T[tuple(i)][1]) if c.converged(tuple(i))
                      else (np.nan, T[tuple(i)][1]))
    start = (LAD[0].index(32), 0, LAD[2].index(1e-6), 0)
    cfg = dict(V121)
    r = run_protocol(look, LAD, start, eps, hs=c.HS, kinds=c.KINDS, **cfg)
    path = list(r["trace"])
    final = tuple(r["proposal"])
    cost_pass = final != tuple(path[-1])
    if cost_pass:
        path.append(final)
    steps = []
    for n, a in enumerate(path):
        lb = LocalBudget(look, LAD, a, hs=c.HS, kinds=c.KINDS, local_order=True,
                         resid_margin=cfg["resid_margin"], resolution=cfg["res_frac"] * eps,
                         cancel_rule=cfg["cancel_rule"], joint_res=cfg["joint_res"],
                         fs_cap=1.25, models=c.models)
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        if np.isfinite(lb.E0) and abs(lb.E0) <= eps:
            lb.joint_correct(look, eps, cancel_rule=cfg["cancel_rule"])
        steps.append(dict(
            setting=[LAD[k][a[k]] for k in range(4)],
            budget=[abs(x) / eps if np.isfinite(x) else None for x in lb.rk],
            estimate=abs(lb.E0) / eps if np.isfinite(lb.E0) else None,
            true=abs(T[a][0][ob] - c.exact[ob]) / eps))
    est = [s["estimate"] for s in steps]
    p1_end = next((n for n, v in enumerate(est) if v is not None and v <= 1.0), len(steps) - 1)
    return dict(target_rel=3e-2, accepted=bool(r["certified"]), runs=int(r["runs"]),
                cost_pass=cost_pass, phase1_end=p1_end, steps=steps)


def l5():
    P = grids.L5
    t = grids.load(P)
    ex = grids.exact(P)["y_min"]
    N = P["lad"][0]
    fin = (len(P["lad"][1]) - 1, len(P["lad"][2]) - 1, len(P["lad"][3]) - 1)
    e = [abs(t[(a,) + fin][0]["y_min"] - ex) / abs(ex) for a in range(len(N))]
    return dict(N=N, err=e, target=3e-4)


if __name__ == "__main__":
    d = dict(dfg=dfg_mesh(), search=search(), l5=l5())
    json.dump(d, open(os.path.join(HERE, "figures", "data_2d1_l5.json"), "w"), indent=1)
    print("written", os.path.join(HERE, "figures", "data_2d1_l5.json"))
