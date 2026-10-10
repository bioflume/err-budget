"""If it is not the slack, what fails on a closely spaced ladder?

Section 7.5 found reliability falling from 100% to 80% as the rungs close from
2x to 1.2x, but the slack -- how far inside the target the recommendation lands
-- is about 0.6 on both, so the first explanation offered was wrong.  The
alternative is that the ESTIMATE degrades: the offset power-law fit of (4) sees
h varying by 16x across a coarse ladder and by only 2x across a fine one, so it
is far less determined there.  This measures both.
"""
import itertools, os, pickle, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from knobs import LocalBudget
from knobs.problems import Burgers1D

D = os.path.join(os.path.dirname(__file__), "data")
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
for rd, rt in [(2.0, 100.0), (1.6, 10.0), (1.4, 5.0), (1.2, 3.0)]:
    path = os.path.join(D, f"burgers_sp_{rd}_{rt}.pkl")
    if not os.path.exists(path):
        print(f"  {rd}x: tensor not built yet"); continue
    tensor = pickle.load(open(path, "rb"))
    LAD = ladders(rd, rt)
    HS = [[2*np.pi/n for n in LAD[0]], [1.0/s for s in LAD[1]],
          list(LAD[2]), list(LAD[3])]
    hr = HS[0][0] / HS[0][-1]
    for ob in ["energy", "dissipation"]:
        QS = EX[ob]
        look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
        rat, fq = [], []
        for a in ALL:
            Et = abs(tensor[a][0][ob] - QS)
            if Et <= abs(QS) * 1e-14: continue
            lb = LocalBudget(look, LAD, a, hs=HS, fs_cap=1.25)
            lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
            if not np.isfinite(lb.E0): continue
            rat.append(abs(lb.E0) / Et)
            good = [q for q in lb.fit_q if np.isfinite(q)]
            if good: fq.append(float(np.median(good)))
        r = np.array(rat)
        print(f"  {f'{rd}x / {rt}x':<18} {hr:>8.1f}x {ob:<13} "
              f"{np.percentile(r,10):>6.2f} /{np.median(r):>6.2f} /"
              f"{np.percentile(r,90):>7.2f} {100*np.mean(r<1):>10.1f}% "
              f"{np.median(fq) if fq else float('nan'):>10.3f}")
print()
print("  'optimistic' is the fraction of anchors where the margined estimate is")
print("  smaller than the true error.  'fit resid' is the median normalised")
print("  residual of the offset power-law fit, the method's own diagnostic.")
