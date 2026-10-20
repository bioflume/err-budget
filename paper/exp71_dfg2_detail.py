"""DFG 2D-2: every target, including those below ten times the reference
uncertainty (not admitted to the pooled table), for the configurations that
differ in how interacting errors are treated (referee points 4 and 8).
    python exp71_dfg2_detail.py

Version 7 of the article: the released program with cancellation credited
(v1.2.1-credit) is added, and for the drag at 3e-3 the script prints the budget
at the setting that this configuration returns outside the target (Section 7.2
and the Remark of Section 4.8).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import cases as C
import pooled as P
from knobs import LocalBudget

case = C._dfg2()
print("reference uncertainty (relative):", {k: f"{v:.1e}" for k, v in case.unc.items()})
for ob, tgs in case.targets.items():
    for tg in tgs:
        adm = "admitted" if case.admitted(ob, tg) else "NOT admitted"
        for cfg in ("v1.1", "v1.1+gated", "v1.2-credit", "v1.2-nojoint", "v1.2", "v1.2.1-credit", "v1.2.1"):
            s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
            print(f"{ob:<7} {tg:<7.0e} {adm:<13} {cfg:<13} starts {s['n']} start-failed {s['start_failed']} "
                  f"accepted {s['accepted']} missed {s['missed']} declined(reach) {s['declined_reach']} "
                  f"pen {s['pen_med']:.2f} worst {s['worst']:.2f}", flush=True)

# The drag at 3e-3 with cancellation credited: the settings returned outside the target and
# the budget there, with and without credit for cancellation (version 1.2.1 otherwise).
ob, tg = "cd_max", 3e-3
o = P.replay(case, ob, tg, P.CONFIGS["v1.2.1-credit"])
eps = tg * abs(case.exact[ob])
T = case.tensor


def look(i):
    v = T[tuple(i)]
    if not np.isfinite(v[1]) or not case.converged(tuple(i)):
        return np.nan, (v[1] if np.isfinite(v[1]) else np.inf)
    return v[0][ob], v[1]


print(f"\n{ob} at {tg:g}, cancellation credited: settings returned outside the target "
      f"(setting: number of starts) {o.get('proposals_missed', {})}")
for q in sorted(o.get("proposals_missed", {})):
    for rule in ("credit", "none"):
        lb = LocalBudget(look, case.LAD, q, hs=case.HS, fs_cap=1.25, resid_margin=5.0, models=case.models,
                         kinds=case.KINDS, local_order=True, resolution=0.1 * eps, cancel_rule=rule,
                         joint_res=1.0)
        lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
        before = lb.E0
        lb.joint_correct(look, eps, cancel_rule=rule)
        print(f"   {[case.LAD[k][j] for k, j in enumerate(q)]} rule={rule}: true error/eps "
              f"{(T[q][0][ob] - case.exact[ob]) / eps:+.3f}, r_k/eps {np.round(np.array(lb.rk) / eps, 3).tolist()}, "
              f"estimate/eps {before / eps:+.3f} before and {lb.E0 / eps:+.3f} after the joint check")
