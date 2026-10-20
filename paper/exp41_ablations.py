"""Table 4: what each design choice buys, measured by switching it off.

All fourteen attainable cells of L2, L3 and L4.  Every
start on each grid is replayed with the method and with one design choice
replaced by its alternative.  The screen and the audit only act when screening
is on, so those two comparisons are made with screening enabled on both sides.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from errbudget import run_protocol
from grids import L2, L3, L4, hs, load, exact, grid

VARIANTS = [
    ("method", dict(use_screen=False)),
    ("propose predicted optimum", dict(use_screen=False, walk_mode="propose")),
    ("Phase 2 coarsen only", dict(use_screen=False, phase2="coarsen")),
    ("method, screened", dict(use_screen=True)),
    ("screen at user's setting", dict(use_screen=True, screen_at="user")),
    ("audit uses one difference", dict(use_screen=True, audit_factor=1.0)),
]

rows = {v[0]: dict(cert=0, safe=0, n=0, pens=[], cap=0, cyc=0, worst_over=0.0)
        for v in VARIANTS}
for P in (L2, L3, L4):
    t, ex, H, ALL = load(P), exact(P), hs(P), grid(P)
    for ob in P["obs"]:
        look = lambda i, ob=ob: (t[tuple(i)][0][ob], t[tuple(i)][1])
        for tg in P["targets"]:
            eps = tg * abs(ex[ob]) if P["rel"] else tg
            err = {i: abs(t[i][0][ob] - ex[ob]) for i in ALL}
            copt = min(t[i][1] for i in ALL if err[i] <= eps)
            for name, kw in VARIANTS:
                R = rows[name]
                for a in ALL:
                    r = run_protocol(look, P["lad"], a, eps, hs=H, kinds=P["kinds"], **kw)
                    R["n"] += 1; R["cap"] += r["hit_cap"]; R["cyc"] += r["cycled"]
                    if r["certified"]:
                        q = r["proposal"]; R["cert"] += 1
                        R["safe"] += err[q] <= eps
                        R["pens"].append(t[q][1] / copt)
                        R["worst_over"] = max(R["worst_over"], err[q] / eps)
            print(f"  done {P['name']} {ob} {tg:g}", flush=True)

print("\n" + "=" * 100)
print(f"{'variant':<28}{'replays':>8}{'certified':>10}{'safety':>9}{'missed':>8}"
      f"{'med pen':>9}{'worst pen':>10}{'hit cap':>9}{'cycled':>8}{'worst err/eps':>14}")
for name, _ in VARIANTS:
    R = rows[name]
    p = np.array(R["pens"]) if R["pens"] else np.array([np.nan])
    print(f"{name:<28}{R['n']:>8}{R['cert']:>10}{100*R['safe']/max(R['cert'],1):>8.1f}%"
          f"{R['cert']-R['safe']:>8}{np.median(p):>9.2f}{p.max():>10.2f}"
          f"{R['cap']:>9}{R['cyc']:>8}{R['worst_over']:>14.2f}")
