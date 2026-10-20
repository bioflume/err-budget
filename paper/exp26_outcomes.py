"""The complete outcome accounting demanded by the referee.

Replaces the summary of exp13 with the quantities a reader needs in order to
judge the reliability claim rather than take it:

  * the full 2x2 outcome table, {certified, refused} x {attainable, not},
    including targets deliberately set below anything the grid can reach, so
    that the false-refusal rate has a denominator;
  * exact replay and failure counts, not a percentage;
  * the number of DISTINCT recommended configurations, which is the honest
    measure of how much independent information an exhaustive tensor replay
    contains;
  * worst case beside the median, since the worst case is what a user risks;
  * the robust-optimum penalty for EVERY cell, not only the cell where the
    ordinary optimum is a cancellation coincidence;
  * safety restricted to the two-point-fit corner -- starts at the coarse end
    of a ladder, where the offset power law is determined exactly and carries
    no residual.  That is where real users start.
"""
import itertools
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from errbudget import run_protocol
from errbudget.problems import Burgers1D, Burgers2D, ReactionDiffusion2D

D = os.path.join(os.path.dirname(__file__), "data")


def neighbours_finer(i, LAD):
    """The one-step-finer neighbours of i, one per ladder."""
    out = []
    for k in range(len(LAD)):
        if i[k] + 1 < len(LAD[k]):
            out.append(tuple(i[k] + 1 if j == k else i[j] for j in range(len(LAD))))
    return out


def cell(tensor, LAD, HS, KIN, ob, QS, t, rel, fs=1.25):
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    eps = t * abs(QS) if rel else t
    err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
    cost = {i: tensor[i][1] for i in ALL}
    feas = [i for i in ALL if err[i] <= eps]
    attainable = bool(feas)

    # robust: attains, and so does every one-step-finer neighbour.  This is the
    # check a careful user performs, and it is the standard the method is held
    # to everywhere rather than only where it flatters.
    robust = [i for i in feas if all(err[j] <= eps for j in neighbours_finer(i, LAD))]
    c_opt = min((cost[i] for i in feas), default=np.nan)
    c_rob = min((cost[i] for i in robust), default=np.nan)

    lk = lambda i: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    n_cert = n_ref = n_safe = n_unsafe = 0
    pens, pens_r, saved, props = [], [], [], []
    coarse_tot = coarse_safe = coarse_cert = 0
    fixed = nbad = 0
    for a in ALL:
        r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KIN, use_screen=False,
                         fs_cap=fs)
        q = r["proposal"]
        at_coarse_end = any(a[k] == 0 for k in range(len(LAD)))
        coarse_tot += at_coarse_end
        if not r["certified"]:
            n_ref += 1
            coarse_cert += 0
            continue
        n_cert += 1
        ok = err[q] <= eps
        n_safe += ok
        n_unsafe += (not ok)
        props.append(q)
        if attainable:
            pens.append(cost[q] / c_opt)
            if robust:
                pens_r.append(cost[q] / c_rob)
        if err[a] <= eps:
            saved.append(cost[a] / cost[q])
        else:
            nbad += 1
            fixed += ok
        if at_coarse_end:
            coarse_cert += 1
            coarse_safe += ok

    return dict(
        observable=ob, target=t, relative=rel, attainable=attainable,
        n_starts=len(ALL), n_certified=n_cert, n_refused=n_ref,
        n_safe=n_safe, n_unsafe=n_unsafe,
        n_attaining=len(feas), n_robust=len(robust),
        distinct_proposals=len(set(props)),
        pen_med=float(np.median(pens)) if pens else np.nan,
        pen_max=float(np.max(pens)) if pens else np.nan,
        pen_rob_med=float(np.median(pens_r)) if pens_r else np.nan,
        pen_rob_max=float(np.max(pens_r)) if pens_r else np.nan,
        saved_med=float(np.median(saved)) if saved else np.nan,
        saved_min=float(np.min(saved)) if saved else np.nan,
        n_already_ok=len(saved), n_missed=nbad, n_fixed=fixed,
        coarse_starts=coarse_tot, coarse_certified=coarse_cert,
        coarse_safe=coarse_safe)


def study(name, tensor, LAD, HS, KIN, exact, obs_list, targets, rel):
    print(f"\n{'='*104}\n{name}\n{'='*104}")
    print(f"{'observable':<12} {'target':<9} {'cert':>9} {'refused':>8} "
          f"{'unsafe':>7} {'distinct':>9} {'pen med/max':>13} "
          f"{'robust med/max':>15} {'attain/robust':>14}")
    rows = []
    for ob in obs_list:
        for t in targets:
            r = cell(tensor, LAD, HS, KIN, ob, exact[ob], t, rel)
            rows.append(dict(problem=name.split()[0], **r))
            tag = "" if r["attainable"] else "  (UNATTAINABLE)"
            print(f"{ob:<12} {t:<9.0e} {r['n_certified']:>5}/{r['n_starts']:<3} "
                  f"{r['n_refused']:>8} {r['n_unsafe']:>7} "
                  f"{r['distinct_proposals']:>9} "
                  f"{r['pen_med']:>6.2f}/{r['pen_max']:<6.2f} "
                  f"{r['pen_rob_med']:>7.2f}/{r['pen_rob_max']:<7.2f} "
                  f"{r['n_robust']:>6}/{r['n_attaining']:<6}{tag}")
    return rows


rows = []

# ------------------------------------------------------------------ L2 ------
b1 = Burgers1D(nu=0.05, T=1.0)
t1 = pickle.load(open(os.path.join(D, "burgers_tensor.pkl"), "rb"))
LAD1 = [[32, 64, 128, 256, 512], [25, 50, 100, 200, 400],
        [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
HS1 = [[2 * np.pi / n for n in LAD1[0]], [1.0 / s for s in LAD1[1]],
       list(LAD1[2]), list(LAD1[3])]
KIN1 = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
rows += study("L2 1D viscous Burgers (Cole-Hopf truth), 4 parameters, 625 starts",
              t1, LAD1, HS1, KIN1, b1.exact_observables(),
              ["energy", "dissipation"], [1e-2, 3e-3, 1e-3, 1e-5], False)

# ------------------------------------------------------------------ L3 ------
b2 = Burgers2D(nu=0.05, T=1.0)
t2 = pickle.load(open(os.path.join(D, "burgers2d_tensor.pkl"), "rb"))
LAD2 = [[16, 32, 64, 128], [10, 20, 40, 80], [1e-2, 1e-5, 1e-8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS2 = [[2 * np.pi / n for n in LAD2[0]], [1.0 / s for s in LAD2[1]],
       list(LAD2[2]), list(LAD2[3]), list(LAD2[4])]
KIN2 = ["discretisation", "discretisation", "nested-solver", "nested-solver",
        "cost-only"]
rows += study("L3 2D Burgers (manufactured solution), 5 parameters, 432 starts",
              t2, LAD2, HS2, KIN2, b2.exact_observables(),
              ["energy", "gradient"], [1e-2, 3e-3, 1e-6], True)

# ------------------------------------------------------------------ L4 ------
p4 = ReactionDiffusion2D()
t4 = pickle.load(open(os.path.join(D, "reactdiff_tensor.pkl"), "rb"))
LAD4 = [[16, 32, 64, 128], [5, 10, 20, 40], [1, 2, 4, 8],
        [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
HS4 = [[2 * np.pi / n for n in LAD4[0]], [1.0 / s for s in LAD4[1]],
       [1.0 / r for r in LAD4[2]], list(LAD4[3]), list(LAD4[4])]
KIN4 = ["discretisation", "discretisation", "discretisation",
        "nested-solver", "cost-only"]
rows += study("L4 Gray-Scott, Strang splitting, 5 parameters, 576 starts",
              t4, LAD4, HS4, KIN4, p4.exact_observables(),
              ["inhibitor_mass", "pattern_sharpness"], [1e-2, 3e-3, 1e-7], True)

json.dump(rows, open(os.path.join(D, "outcomes.json"), "w"), indent=1)

# ----------------------------------------------------------------- totals ---
att = [r for r in rows if r["attainable"]]
un = [r for r in rows if not r["attainable"]]
print(f"\n{'='*104}\nTOTALS\n{'='*104}")
print(f"  cells with the target attainable on the grid : {len(att)}")
print(f"  cells with the target below anything the grid can reach : {len(un)}")
print()
print("                                        certified      refused")
print(f"  target attainable on the grid      {sum(r['n_certified'] for r in att):9d} "
      f"{sum(r['n_refused'] for r in att):12d}")
print(f"  target NOT attainable              {sum(r['n_certified'] for r in un):9d} "
      f"{sum(r['n_refused'] for r in un):12d}")
print()
ncert = sum(r["n_certified"] for r in rows)
nun = sum(r["n_unsafe"] for r in rows)
nrep = sum(r["n_starts"] for r in rows)
print(f"  replays graded                     : {nrep}")
print(f"  certified                          : {ncert}")
print(f"  certified and target NOT attained  : {nun}")
print(f"  false refusals (attainable cell)   : {sum(r['n_refused'] for r in att)}"
      f"  = {100*sum(r['n_refused'] for r in att)/max(sum(r['n_starts'] for r in att),1):.1f}%"
      f" of starts in attainable cells")
if nun == 0:
    print(f"  one-sided 95% upper bound on the failure rate, treating each replay")
    print(f"  as independent (it is not): 3/{ncert} = {300/ncert:.3f}%")
    print(f"  the same bound at the honest unit of independence, the cell: "
          f"3/{len(att)} = {300/len(att):.0f}%")
print()
cs = sum(r["coarse_safe"] for r in att)
cc = sum(r["coarse_certified"] for r in att)
print(f"  two-point-fit corner (a start at the coarse end of some ladder):")
print(f"    {cc} certified, {cs} attained the target, "
      f"{cc-cs} did not")
