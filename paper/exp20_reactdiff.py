"""L4: does the method work on a different domain and a different architecture?

Gray-Scott reaction-diffusion, Strang splitting, no Newton iteration, three
competing discretisation errors.  Replayed from every starting configuration.
"""
import sys, os, itertools, pickle, time, numpy as np
from knobs import run_protocol
from knobs.problems import ReactionDiffusion2D

NAMES = ["N", "steps", "n_react", "gmres_tol", "ilu_drop"]
LAD   = [[16, 32, 64, 128], [5, 10, 20, 40], [1, 2, 4, 8],
         [1e-2, 1e-5, 1e-8], [1e-1, 1e-2, 1e-3]]
KINDS = ["discretisation", "discretisation", "discretisation",
         "nested-solver", "cost-only"]
HS = [[2*np.pi/n for n in LAD[0]], [1.0/s for s in LAD[1]],
      [1.0/r for r in LAD[2]], list(LAD[3]), list(LAD[4])]
ALL = list(itertools.product(*[range(len(l)) for l in LAD]))

p = ReactionDiffusion2D()
EXACT = p.exact_observables()
CACHE = os.path.join(os.path.dirname(__file__), "data", "reactdiff_tensor.pkl")
if os.path.exists(CACHE):
    tensor = pickle.load(open(CACHE, "rb"))
else:
    tensor, t0 = {}, time.time()
    for i, idx in enumerate(ALL):
        th = {NAMES[j]: LAD[j][v] for j, v in enumerate(idx)}
        obs, C, w, _ = p.run(th)
        tensor[idx] = (obs, C)
        if i % 100 == 0:
            print(f"  {i}/{len(ALL)}  {time.time()-t0:.0f}s", flush=True)
    pickle.dump(tensor, open(CACHE, "wb"))
    print(f"tensor: {len(tensor)} runs in {time.time()-t0:.0f}s")

def show(i):
    return "N=%d,steps=%d,n_react=%d,gt=%.0e,dr=%.0e" % tuple(
        LAD[k][j] for k, j in enumerate(i))

print(f"\n{'='*94}\nL4  2D Gray-Scott reaction-diffusion, 5 params, {len(ALL)} starting points")
print("observable          target     meets   penalty(med/p90)  already-ok: saved  missed: fixed  calib")
for ob in ["inhibitor_mass", "pattern_sharpness"]:
    QS = EXACT[ob]
    look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    tr = lambda i: abs(tensor[i][0][ob] - QS)
    for er in [1e-2, 3e-3]:
        eps = er * abs(QS)
        feas = [i for i in ALL if tr(i) <= eps]
        if not feas:
            print(f"{ob:<19} {er:.0e}rel  infeasible (best {min(tr(i) for i in ALL)/abs(QS):.2e} rel)")
            continue
        c_opt = min(tensor[i][1] for i in feas)
        safe = 0; pens = []; saved = []; nbad = 0; fixed = 0; cal = []
        for a in ALL:
            r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, use_screen=False)
            q = r["proposal"]; ok = tr(q) <= eps; safe += ok
            pens.append(tensor[q][1] / c_opt); cal.append(r["protocol_cost"]/tensor[a][1])
            if tr(a) <= eps: saved.append(tensor[a][1] / tensor[q][1])
            else: nbad += 1; fixed += ok
        pens = np.array(pens); sv = np.array(saved) if saved else np.array([np.nan])
        print(f"{ob:<19} {er:.0e}rel  {100*safe/len(ALL):5.1f}%  "
              f"{np.median(pens):5.2f}/{np.percentile(pens,90):5.2f}     "
              f"{np.nanmedian(sv):6.2f}x ({len(sv):3d})    {fixed:3d}/{nbad:<3d}     "
              f"{np.median(cal):6.1f}x")

print(f"\n{'='*94}\nerror budget at a mid configuration (what the user is told)")
from knobs import LocalBudget
for ob in ["inhibitor_mass", "pattern_sharpness"]:
    look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
    lb = LocalBudget(look, LAD, (2, 2, 2, 1, 1), hs=HS)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    tot = sum(abs(r) for r in lb.rk if np.isfinite(r)) or 1.0
    print(f"  {ob}:")
    for k, nm in enumerate(NAMES):
        print(f"    {nm:<12} order {lb.pk[k]:4.2f}   residual {lb.rk[k]:+.3e}   "
              f"{100*abs(lb.rk[k])/tot:5.1f}% of the budget")
