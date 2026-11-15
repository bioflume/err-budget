"""Check for the spacing study (exp34_spacing.py): at spacing 1.2 / 3, does the
misfit-dependent safety factor cause the declines?  Replays the released program
(errbudget.V121) with and without it (resid_margin=0, a constant factor 1.25) on the
same four test cases, graded as in exp34_spacing.py.  Results go to
data/spacing_margin_check.json.
"""
import itertools
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from errbudget import run_protocol, V121
from errbudget.problems import Burgers1D
from l2flags import flags
from paths import DATA, OUT

KIN = ["discretisation", "discretisation", "nested-solver", "nested-solver"]
ALL = list(itertools.product(*[range(5) for _ in range(4)]))
EX = Burgers1D(nu=0.05, T=1.0).exact_observables()
rd, rt = 1.2, 3.0
LAD = [[int(round(128 * rd ** (j - 2))) for j in range(5)],
       [int(round(100 * rd ** (j - 2))) for j in range(5)],
       [1e-6 * rt ** (2 - j) for j in range(5)], [1e-6 * rt ** (2 - j) for j in range(5)]]
HS = [[2 * np.pi / n for n in LAD[0]], [1.0 / s for s in LAD[1]], list(LAD[2]), list(LAD[3])]
tensor = pickle.load(open(os.path.join(DATA, f"burgers_sp_{rd}_{rt}.pkl"), "rb"))
FL = flags(tensor, LAD, 0.05, f"L2_sp_{rd}_{rt}")
out = {}
for name, kw in [("released", dict(V121)), ("released, constant factor", dict(V121, resid_margin=0.0))]:
    acc = met = 0
    for ob in ["energy", "dissipation"]:
        QS = EX[ob]
        err = {i: abs(tensor[i][0][ob] - QS) for i in ALL}
        lk = lambda i, ob=ob: ((tensor[tuple(i)][0][ob] if FL[tuple(i)]["converged"] else np.nan),
                               tensor[tuple(i)][1])
        cand = sorted(err.values())
        for eps in [cand[int(0.25 * len(cand))], cand[int(0.5 * len(cand))]]:
            for a in ALL:
                r = run_protocol(lk, LAD, a, eps, hs=HS, kinds=KIN, **kw)
                if r["certified"]:
                    acc += 1
                    met += err[r["proposal"]] <= eps
    out[name] = dict(accepted=acc, met=met, searches=4 * len(ALL))
    print(f"spacing {rd}/{rt}, {name}: accepted {acc}, met {met} of {4 * len(ALL)}", flush=True)
json.dump(out, open(os.path.join(OUT, "spacing_margin_check.json"), "w"), indent=1)
