"""Analysis of J1: function object against the article's post-processing, drag
convergence to n_arc = 192 with a current release, and the program's results on
the 2D-1 grid computed with the current release.
    python analyse_j1.py [config ...]
"""
import json
import os
import pickle
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "paper"))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import dfg2d1_hpc as J
from openfoam_cylinder import REF
import cases as C
import pooled as P

RES = os.path.join(os.environ.get("WORK", "."), "dfg2d1_hpc")
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

print("forces function object against the article's post-processing (current release)")
print(f"{'n_arc':>6} {'iters':>6} {'Cd (FO)':>12} {'Cd (fields)':>12} {'Cd (v1.1)':>12}"
      f" {'err FO':>9} {'err fields':>10} {'err v1.1':>9} {'Cl FO':>10} {'Cl fields':>10}")
fo = {}
for n in J.FO_LEVELS:
    f = os.path.join(RES, f"fo_na{n}", "result.json")
    if not os.path.exists(f):
        continue
    r = json.load(open(f)); fo[n] = r
    print(f"{n:>6} {r['iters']:>6} {r.get('cd_fo', np.nan):12.6f} {r['cd_fields']:12.6f} "
          f"{r['cd_fields_vertexmean']:12.6f} {r.get('err_cd_fo', np.nan):9.2e} "
          f"{r['err_cd_fields']:10.2e} {r['err_cd_fields_vertexmean']:9.2e} "
          f"{r.get('cl_fo', np.nan):10.6f} {r['cl_fields']:10.6f}")
json.dump(fo, open(os.path.join(OUT, "j1_fo.json"), "w"), indent=1)

tensor = {}
for t in J.TASKS():
    if t["kind"] != "grid":
        continue
    f = os.path.join(RES, t["id"], "result.json")
    if os.path.exists(f):
        r = json.load(open(f))
        tensor[tuple(r["idx"])] = (r["obs"], r["cost"], r["info"])
print(f"\ngrid with the current release: {len(tensor)}/840 settings")
if len(tensor) == 840:
    pickle.dump(tensor, open(os.path.join(OUT, "dfg_tensor_current.pkl"), "wb"))
    base = C._dfg1()
    cur = C.Case("2D-1c", "OpenFOAM DFG 2D-1, current release", {}, base.LAD, base.HS,
                 base.KINDS, base.exact, base.targets)
    for i, v in tensor.items():
        info = dict(v[2]); c = v[1]
        if np.isfinite(c) and info.get("cells") and info.get("iters"):
            c = float(info["cells"] * info["iters"])
        cur.tensor[i] = (v[0], c, info)
    rows = []
    for cfg in (sys.argv[1:] or ["v1.2"]):
        for case in (base, cur):
            for ob, tgs in case.targets.items():
                for tg in tgs:
                    s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
                    s.update(config=cfg, release="v1912" if case is base else "current")
                    rows.append(s)
                    print(f"  {cfg} {s['release']:<8} {ob:<3} {tg:.0e}: accepted {s['accepted']}"
                          f"/{s['n']} missed {s['missed']} penalty {s['pen_med']:.2f}", flush=True)
    json.dump(rows, open(os.path.join(OUT, "j1_replay.json"), "w"), indent=1)
