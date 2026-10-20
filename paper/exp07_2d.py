"""L3: 2D Burgers, five knobs, one of them accuracy-neutral."""
import sys, os, itertools, pickle, time, os, numpy as np
from errbudget import LocalBudget, run_protocol
from errbudget.problems import Burgers2D

b = Burgers2D(nu=0.05, T=1.0)
EXACT = b.exact_observables()
OBS = ["energy", "gradient"]

Ns = [16, 32, 64, 128]
Ss = [10, 20, 40, 80]
NT = [1e-2, 1e-5, 1e-8]
GT = [1e-2, 1e-5, 1e-8]
DR = [1e-1, 1e-2, 1e-3]
LAD = [Ns, Ss, NT, GT, DR]
NAMES = ["N", "steps", "newton_tol", "gmres_tol", "ilu_drop"]
HS = [[2 * np.pi / n for n in Ns], [1.0 / s for s in Ss], list(NT), list(GT), list(DR)]
ALL = list(itertools.product(*[range(len(l)) for l in LAD]))

CACHE = os.path.join(os.path.dirname(__file__),"data","burgers2d_tensor.pkl")
if os.path.exists(CACHE):
    tensor = pickle.load(open(CACHE, "rb"))
else:
    tensor, t0 = {}, time.time()
    for n, idx in enumerate(ALL):
        th = dict(zip(NAMES, [LAD[i][j] for i, j in enumerate(idx)]))
        obs, C, w, _ = b.run(th)
        tensor[idx] = (obs, C, w)
        if n % 50 == 0:
            print(f"  {n}/{len(ALL)}  {time.time() - t0:.0f}s", flush=True)
    pickle.dump(tensor, open(CACHE, "wb"))
    print(f"tensor: {len(tensor)} runs in {time.time() - t0:.0f}s")


def show(idx):
    return "N=%d,st=%d,nt=%.0e,gt=%.0e,dr=%.0e" % tuple(
        LAD[i][j] for i, j in enumerate(idx))


for ob in OBS:
    QS = EXACT[ob]
    look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    tr = lambda i: abs(tensor[i][0][ob] - QS)
    print(f"\n{'=' * 86}\n=== {ob}  (exact {QS:.8f}) ===")

    # what does the sweep say about each knob, at a mid anchor?
    mid = (2, 2, 1, 1, 1)
    lb = LocalBudget(look, LAD, mid, hs=HS)
    lb.sweep(); lb.per_knob_error(); lb.refine_path()
    print(f"  at {show(mid)}:  estimated error {lb.E0:+.3e}  (true {tensor[mid][0][ob]-QS:+.3e})")
    print("    knob          p_hat   residual r_k    budget share   masked?")
    thr = lb.thresholds()
    tot = sum(abs(r) for r in lb.rk if np.isfinite(r)) or 1.0
    for k, nm in enumerate(NAMES):
        print(f"    {nm:<12} {lb.pk[k]:6.2f}   {lb.rk[k]:+.3e}    {100*abs(lb.rk[k])/tot:5.1f}%"
              f"         {thr[k]['masked']}")

    for eps_rel in [1e-2, 3e-3]:
        eps = eps_rel * abs(QS)
        feas = [i for i in ALL if tr(i) <= eps]
        if not feas:
            print(f"\n  rel target {eps_rel:.0e}: infeasible")
            continue
        opt = min(feas, key=lambda i: tensor[i][1])
        c_opt = tensor[opt][1]
        # A user who was already accurate wants a cheaper run; a user who was
        # not wants a correct one.  Averaging the two produces a "saving" that
        # describes neither, so split them.
        safe, pens, ovh, runs = 0, [], [], []
        save_ok, cost_up_bad, n_bad_fixed, n_bad = [], [], 0, 0
        for a in ALL:
            r = run_protocol(look, LAD, a, eps, hs=HS)
            p = r["proposal"]
            ok = tr(p) <= eps
            safe += ok
            pens.append(tensor[p][1] / c_opt)
            ovh.append(r["protocol_cost"] / tensor[a][1])
            runs.append(r["runs"])
            if tr(a) <= eps:                       # user was already accurate
                save_ok.append(tensor[a][1] / tensor[p][1])
            else:                                  # user's run missed the target
                n_bad += 1
                n_bad_fixed += ok
                cost_up_bad.append(tensor[p][1] / tensor[a][1])
        pens = np.array(pens)
        print(f"\n  rel target {eps_rel:.0e} (abs {eps:.3e})   true optimum {show(opt)}")
        print(f"    safe ................. {safe}/{len(ALL)} ({100*safe/len(ALL):.1f}%)")
        print(f"    penalty vs optimum ... median {np.median(pens):.2f}x  "
              f"p90 {np.percentile(pens,90):.2f}x  max {pens.max():.2f}x")
        if save_ok:
            s = np.array(save_ok)
            print(f"    starts already accurate ({len(s)}): cost saved median "
                  f"{np.median(s):.2f}x  p90 {np.percentile(s,90):.2f}x  max {s.max():.1f}x")
        if n_bad:
            c = np.array(cost_up_bad)
            print(f"    starts that missed target ({n_bad}): fixed {n_bad_fixed}/{n_bad}, "
                  f"cost median {np.median(c):.2f}x more")
        print(f"    calibration .......... median {np.median(ovh):.1f}x the user's own run, "
              f"{np.median(runs):.0f} distinct runs")

# ------------------------------------------------- the accuracy-neutral knob
print(f"\n{'=' * 86}\n=== does the protocol find the cost minimum of the neutral knob? ===")
for ob in OBS:
    QS = EXACT[ob]
    look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    tr = lambda i: abs(tensor[i][0][ob] - QS)
    eps = 1e-2 * abs(QS)
    hit, chosen, best_idx = 0, [], []
    for a in ALL:
        r = run_protocol(look, LAD, a, eps, hs=HS)
        p = r["proposal"]
        # the right baseline is the cheapest drop tolerance holding THIS
        # proposal's other four knobs fixed, not the global optimum's
        costs = [tensor[tuple(p[j] if j != 4 else d for j in range(5))][1]
                 for d in range(len(DR))]
        bi = int(np.argmin(costs))
        hit += (p[4] == bi)
        chosen.append(p[4]); best_idx.append(bi)
        if a == ALL[0]:
            excess = costs[p[4]] / costs[bi]
    chosen, best_idx = np.array(chosen), np.array(best_idx)
    exc = []
    for a, c, bi in zip(ALL, chosen, best_idx):
        pass
    print(f"  {ob:<10} picks the cost-minimising drop tolerance in {hit}/{len(ALL)} "
          f"({100*hit/len(ALL):.0f}%).  chosen {np.bincount(chosen, minlength=len(DR))}  "
          f"optimal {np.bincount(best_idx, minlength=len(DR))}  (0=1e-1,1=1e-2,2=1e-3)")
