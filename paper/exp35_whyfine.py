"""If it is not the slack, what fails on a closely spaced ladder?

Section 7.5 found reliability falling from 100% to 80% as the rungs close from
2x to 1.2x, but the slack -- how far inside the target the recommendation lands
-- is about 0.6 on both, so the first explanation offered was wrong.  The
alternative is that the ESTIMATE degrades: the offset power-law fit of (4) sees
h varying by 16x across a coarse ladder and by only 2x across a fine one, so it
is far less determined there.  This measures both.

Rerun for the revision without version history: the estimate is that of the released
program (errbudget.V121, built as in exp66_estimate_v12.py, including the interaction
check), with the resolution set from the tighter of the two targets of exp34_spacing.py;
unconverged runs are failed (flags from l2flags.py).  Reported per spacing and observable:
the median misfit of the fits, and the 10th percentile, median and 90th percentile of
Ehat/Etrue with the share of optimistic estimates, over every setting whose true error
exceeds the resolution (as exp66 grades it).  Results go to data/spacing_estimate_released.json.
"""
import itertools, os, pickle
import numpy as np
import json
from errbudget import LocalBudget, V121
from l2flags import flags
from errbudget.problems import Burgers1D

from paths import DATA as D, OUT
ALL = list(itertools.product(*[range(5) for _ in range(4)]))
p = Burgers1D(nu=0.05, T=1.0)
EX = p.exact_observables()

def ladders(rd, rt):
    N = [int(round(128 * rd ** (j - 2))) for j in range(5)]
    S = [int(round(100 * rd ** (j - 2))) for j in range(5)]
    T = [1e-6 * rt ** (2 - j) for j in range(5)]
    return [N, S, list(T), list(T)]

print("=" * 96)
print("what degrades as the ladder closes up?")
print("=" * 96)
print(f"  {'spacing':<18} {'h range':>9} {'observable':<13} {'Ehat/Etrue':>22} "
      f"{'optimistic':>11} {'fit resid':>10}")
print(f"  {'':<18} {'':>9} {'':<13} {'p10 / median / p90':>22}")
KIN = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
OUTROWS = []
for rd, rt in [(2.0, 100.0), (1.6, 10.0), (1.4, 5.0), (1.2, 3.0)]:
    path = os.path.join(D, f"burgers_sp_{rd}_{rt}.pkl")
    if not os.path.exists(path):
        print(f"  {rd}x: tensor not built yet"); continue
    tensor = pickle.load(open(path, "rb"))
    LAD = ladders(rd, rt)
    FL = flags(tensor, LAD, 0.05, f"L2_sp_{rd}_{rt}")
    HS = [[2*np.pi/n for n in LAD[0]], [1.0/s for s in LAD[1]],
          list(LAD[2]), list(LAD[3])]
    hr = HS[0][0] / HS[0][-1]
    for ob in ["energy", "dissipation"]:
        QS = EX[ob]
        look = lambda i, ob=ob: ((tensor[tuple(i)][0][ob] if FL[tuple(i)]["converged"] else np.nan),
                                 tensor[tuple(i)][1])
        err = sorted(abs(tensor[i][0][ob] - QS) for i in ALL)
        eps = err[int(0.25 * len(err))]          # the tighter target of exp34_spacing.py
        cfg = V121
        rat, fq = [], []
        for a in ALL:
            Et = abs(tensor[a][0][ob] - QS)
            if Et <= 0.1 * eps: continue         # below the resolution, as in exp66
            lb = LocalBudget(look, LAD, a, hs=HS, fs_cap=1.25, resid_margin=cfg["resid_margin"],
                             kinds=KIN, local_order=cfg["local_order_check"],
                             resolution=cfg["res_frac"] * eps, cancel_rule=cfg["cancel_rule"],
                             joint_res=cfg["joint_res"])
            lb.sweep(fine_reach=1); lb.per_knob_error(); lb.refine_path(n_levels=0)
            if np.isfinite(lb.E0):
                lb.joint_correct(look, eps, cancel_rule=cfg["cancel_rule"])
            if not np.isfinite(lb.E0): continue
            rat.append(abs(lb.E0) / Et)
            good = [q for q in lb.fit_q if np.isfinite(q)]
            if good: fq.append(float(np.median(good)))
        r = np.array(rat)
        OUTROWS.append(dict(rd=rd, rt=rt, h_range=hr, ob=ob, eps_rel=eps / abs(QS), n=len(r),
                            p10=float(np.percentile(r, 10)), median=float(np.median(r)),
                            p90=float(np.percentile(r, 90)), optimistic=float(np.mean(r < 1)),
                            misfit=float(np.median(fq)) if fq else None))
        print(f"  {f'{rd}x / {rt}x':<18} {hr:>8.1f}x {ob:<13} "
              f"{np.percentile(r,10):>6.2f} /{np.median(r):>6.2f} /"
              f"{np.percentile(r,90):>7.2f} {100*np.mean(r<1):>10.1f}% "
              f"{np.median(fq) if fq else float('nan'):>10.3f}")
json.dump(OUTROWS, open(os.path.join(OUT, "spacing_estimate_released.json"), "w"), indent=1)
print()
print("  'optimistic' is the fraction of anchors where the margined estimate is")
print("  smaller than the true error.  'fit resid' is the median normalised")
print("  residual of the offset power-law fit, the method's own diagnostic.")
