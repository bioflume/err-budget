"""Is the cost penalty real, or is it the price of not gambling on cancellation?

The "true optimum" found by exhaustive search may rely on error cancellation: it
meets the target, but REFINING it makes the answer worse, so it is accurate by
luck rather than by convergence.  A user who checked their answer the obvious way
-- refine and see whether it moves -- would reject it.  We therefore also grade
against the cheapest ROBUST configuration: one that meets the target and whose
every finer neighbour does too.
"""
import os, itertools, pickle, numpy as np
from errbudget import run_protocol
from errbudget.problems import ReactionDiffusion2D, Burgers2D

def study(name, tensor, LAD, KINDS, HS, obs_list, exact, targets):
    ALL = list(itertools.product(*[range(len(l)) for l in LAD]))
    print(f"\n{'='*92}\n{name}  ({len(ALL)} starting points)")
    print("observable           target   meets   penalty vs best   penalty vs robust   "
          "best relies on luck?")
    for ob in obs_list:
        QS = exact[ob]
        look = lambda i, ob=ob: (tensor[tuple(i)][0][ob], tensor[tuple(i)][1])
        tr = lambda i: abs(tensor[i][0][ob] - QS)
        for er in targets:
            eps = er * abs(QS)
            feas = [i for i in ALL if tr(i) <= eps]
            if not feas:
                continue
            def robust(i):
                for k in range(len(LAD)):
                    if i[k] + 1 < len(LAD[k]):
                        j = tuple(i[k]+1 if m == k else i[m] for m in range(len(LAD)))
                        if tr(j) > eps:
                            return False
                return True
            rob = [i for i in feas if robust(i)]
            c_best = min(tensor[i][1] for i in feas)
            c_rob = min(tensor[i][1] for i in rob) if rob else np.nan
            lucky = "yes" if rob and c_rob > c_best * 1.5 else "no"
            safe, pb, pr = 0, [], []
            for a in ALL:
                r = run_protocol(look, LAD, a, eps, hs=HS, kinds=KINDS, use_screen=False)
                q = r["proposal"]; safe += tr(q) <= eps
                pb.append(tensor[q][1] / c_best)
                if rob:
                    pr.append(tensor[q][1] / c_rob)
            pb = np.array(pb); pr = np.array(pr) if pr else np.array([np.nan])
            print(f"{ob:<20} {er:.0e}   {100*safe/len(ALL):5.1f}%   "
                  f"{np.median(pb):5.2f}x            {np.nanmedian(pr):5.2f}x"
                  f"              {lucky}  ({len(rob)}/{len(feas)} robust)")

rd = ReactionDiffusion2D()
study("L4  2D Gray-Scott reaction-diffusion (splitting, no Newton)",
      pickle.load(open("paper/data/reactdiff_tensor.pkl", "rb")),
      [[16,32,64,128],[5,10,20,40],[1,2,4,8],[1e-2,1e-5,1e-8],[1e-1,1e-2,1e-3]],
      ["discretisation","discretisation","discretisation","nested-solver","cost-only"],
      [[2*np.pi/n for n in [16,32,64,128]],[1/s for s in [5,10,20,40]],
       [1/r for r in [1,2,4,8]],[1e-2,1e-5,1e-8],[1e-1,1e-2,1e-3]],
      ["inhibitor_mass","pattern_sharpness"], rd.exact_observables(), [1e-2, 3e-3])

b2 = Burgers2D(nu=0.05, T=1.0)
study("L3  2D Burgers (Newton + GMRES)",
      pickle.load(open("paper/data/burgers2d_tensor.pkl", "rb")),
      [[16,32,64,128],[10,20,40,80],[1e-2,1e-5,1e-8],[1e-2,1e-5,1e-8],[1e-1,1e-2,1e-3]],
      ["discretisation","discretisation","nested-solver","nested-solver","cost-only"],
      [[2*np.pi/n for n in [16,32,64,128]],[1/s for s in [10,20,40,80]],
       [1e-2,1e-5,1e-8],[1e-2,1e-5,1e-8],[1e-1,1e-2,1e-3]],
      ["energy","gradient"], b2.exact_observables(), [1e-2, 3e-3])
