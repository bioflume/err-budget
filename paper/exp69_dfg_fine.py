"""DFG 2D-1 beyond the grid of the study: n_arc = 96 and 128 at the tightest solver
settings, with the post-processing of version 1.2 (exact wall-cell centroids,
second-order wall pressure), iteration limit 20,000 (referee point 7).
    python exp69_dfg_fine.py
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import openfoam_cylinder as OC
import cases as C

_cd = OC.control_dict
OC.control_dict = lambda end=20000: _cd(end)
OUT = os.path.join(C.D, "dfg_fine.json")
res = json.load(open(OUT)) if os.path.exists(OUT) else {}
p = OC.DFGCylinder()
for na in (64, 96, 128):
    if str(na) in res:
        continue
    obs, cost, wall, info = p.run(dict(n_arc=na, res_ctrl=1e-7, lin_tol=1e-8, n_coarsest=10))
    res[str(na)] = dict(obs=obs, info=info, wall=wall)
    json.dump(res, open(OUT, "w"), indent=1)
    print(na, info, {k: f"{(v - OC.REF[k.split('_')[0]]) / OC.REF[k.split('_')[0]]:+.2e}"
                     for k, v in obs.items()}, flush=True)
