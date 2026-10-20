"""Does the two-anchor protocol work from *any* starting point a user might pick?

We replay the protocol from all 625 configurations in the tensor sweep, so the
result is a distribution over plausible users rather than one anecdote.
"""
import sys, os, itertools, pickle, numpy as np
from errbudget import LocalBudget, run_protocol
from errbudget.problems import Burgers1D

b = Burgers1D(nu=0.05, T=1.0)
EXACT = b.exact_observables()
tensor = pickle.load(open(os.path.join(os.path.dirname(__file__),"data","burgers_tensor.pkl"), "rb"))

Ns, Ss = [32, 64, 128, 256, 512], [25, 50, 100, 200, 400]
NT = GT = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]
LAD = [Ns, Ss, NT, GT]
import numpy as _np
HS = [[2*_np.pi/n for n in Ns], [1.0/s for s in Ss], list(NT), list(GT)]
ALL = list(itertools.product(*[range(5) for _ in range(4)]))


def make_lookup(ob):
    return lambda idx: (tensor[tuple(idx)][0][ob], tensor[tuple(idx)][1])


def truth(ob, idx):
    return abs(tensor[idx][0][ob] - EXACT[ob])


def show(idx):
    return "N=%d,st=%d,nt=%.0e,gt=%.0e" % tuple(LAD[i][j] for i, j in enumerate(idx))


for ob in ["energy", "dissipation"]:
    look = make_lookup(ob)
    QS = EXACT[ob]
    print(f"\n{'=' * 80}\n=== {ob} ===")
    for eps in [1e-2, 3e-3, 1e-3]:
        feas = [i for i in ALL if truth(ob, i) <= eps]
        if not feas:
            print(f"  eps={eps:.0e}: infeasible")
            continue
        c_opt = min(tensor[i][1] for i in feas)
        opt = min(feas, key=lambda i: tensor[i][1])

        safe = 0
        pens, bens, ovh_rec, ovh_anchor, breakeven, iters, unrel = [], [], [], [], [], [], 0
        for a in ALL:
            r = run_protocol(look, LAD, a, eps, mode="optimistic", hs=HS)
            p = r["proposal"]
            safe += truth(ob, p) <= eps
            pens.append(tensor[p][1] / c_opt)
            saved = tensor[a][1] - tensor[p][1]
            bens.append(tensor[a][1] / tensor[p][1])
            ovh_rec.append(r["protocol_cost"] / tensor[p][1])
            ovh_anchor.append(r["protocol_cost"] / tensor[a][1])
            breakeven.append(r["protocol_cost"] / saved if saved > 0 else np.inf)
            iters.append(r["iters"])
            unrel += not r["reliable"]
        pens, bens = np.array(pens), np.array(bens)
        ovh_rec, ovh_anchor = np.array(ovh_rec), np.array(ovh_anchor)
        be = np.array(breakeven)
        print(f"\n  eps={eps:.0e}   true optimum: {show(opt)}  cost {c_opt:.4g}")
        print(f"    safe proposals ............ {safe}/{len(ALL)} "
              f"({100 * safe / len(ALL):.1f}%)")
        print(f"    cost penalty vs optimum ... median {np.median(pens):.2f}x   "
              f"p90 {np.percentile(pens, 90):.2f}x   max {pens.max():.2f}x")
        print(f"    saving vs user's own start  median {np.median(bens):.2f}x   "
              f"p90 {np.percentile(bens, 90):.2f}x   max {bens.max():.1f}x")
        print(f"    calibration cost .......... median {np.median(ovh_rec):.1f}x the "
              f"recommended run, {np.median(ovh_anchor):.2f}x the user's own run")
        fin = be[np.isfinite(be)]
        print(f"    break-even after .......... median {np.median(fin):.1f} production "
              f"runs ({100 * len(fin) / len(be):.0f}% of starts ever break even)")
        print(f"    re-anchoring iterations ... mean {np.mean(iters):.2f}   "
              f"flagged-unreliable: {unrel}/{len(ALL)}")

# ------------------------------------------------- where does it still fail?
print(f"\n{'=' * 80}\n=== failure anatomy (dissipation, eps=1e-3) ===")
ob, eps = "dissipation", 1e-3
look = make_lookup(ob)
bad = []
for a in ALL:
    r = run_protocol(look, LAD, a, eps, mode="optimistic", hs=HS)
    p = r["proposal"]
    if truth(ob, p) > eps:
        bad.append((a, p, truth(ob, p), r["reliable"]))
print(f"  {len(bad)} unsafe proposals out of {len(ALL)}")
if bad:
    print("  starting point                 -> proposal                       "
          "true err   Richardson ok?")
    for a, p, e, rel in bad[:12]:
        print(f"  {show(a):<30} -> {show(p):<30} {e:.2e}   {rel}")
    anch = np.array([a for a, _, _, _ in bad])
    print(f"\n  unsafe starts by knob index (0=coarsest):")
    for i, nm in enumerate(["N", "steps", "newton_tol", "gmres_tol"]):
        print(f"    {nm:<11} {np.bincount(anch[:, i], minlength=5)}")

# ------------------------------------------------ conservative mode as a guard
print(f"\n{'=' * 80}\n=== does the cancellation-free mode buy safety? ===")
for ob in ["energy", "dissipation"]:
    look = make_lookup(ob)
    for eps in [1e-3]:
        feas = [i for i in ALL if truth(ob, i) <= eps]
        c_opt = min(tensor[i][1] for i in feas)
        for mode in ["optimistic", "conservative"]:
            safe, pens = 0, []
            for a in ALL:
                r = run_protocol(look, LAD, a, eps, mode=mode, hs=HS)
                p = r["proposal"]
                safe += truth(ob, p) <= eps
                pens.append(tensor[p][1] / c_opt)
            print(f"  {ob:<12} eps={eps:.0e} {mode:<13} safe {100 * safe / len(ALL):5.1f}%   "
                  f"median penalty {np.median(pens):.2f}x")
