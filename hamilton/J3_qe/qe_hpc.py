"""J3: the Quantum ESPRESSO study with a current release (referee point 26).

The same 480-setting grid, reference and reference uncertainty as the article
(exp53_qe.py), run with the pw.x of the conda environment (QE >= 7.3) instead
of the Ubuntu 24.04 package of pw.x 6.7 (which needed a preloaded shim).
    python -m farm qe_hpc array K NK
    python analyse_j3.py
"""
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
from qe_silicon import QESilicon

NAMES = ["ecutwfc", "nk", "conv_thr", "mixing_beta"]
LAD = [[12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 40.0, 48.0], [2, 3, 4, 6],
       [1e-4, 1e-6, 1e-8, 1e-10, 1e-12], [0.3, 0.5, 0.7]]
SMOKE = os.environ.get("ERRBUDGET_SMOKE") == "1"
_P = {}


def TASKS():
    p0 = dict(ecutwfc=90.0, nk=8, conv_thr=1e-14, mixing_beta=0.7)
    T = [dict(id="ref", theta=p0), dict(id="e75", theta=dict(p0, ecutwfc=75.0)),
         dict(id="k8", theta=dict(p0, ecutwfc=60.0)), dict(id="k10", theta=dict(p0, ecutwfc=60.0, nk=10))]
    for idx in itertools.product(*[range(len(l)) for l in LAD]):
        T.append(dict(id="g_" + "_".join(map(str, idx)), idx=list(idx),
                      theta={NAMES[k]: LAD[k][j] for k, j in enumerate(idx)}))
    for t in T:
        t["nprocs"] = 1
    return T[:6] if SMOKE else T


def run(task, out):
    if "p" not in _P:
        _P["p"] = QESilicon()
    obs, cost, wall, info = _P["p"].run(task["theta"])
    json.dump(dict(task, obs=obs, cost=cost, wall=wall, info=info),
              open(os.path.join(out, "result.json"), "w"), indent=1)
