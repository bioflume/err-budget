"""Accuracy of the error estimate at every setting (Appendix A, Table A.1).

For every usable, converged setting of L2-L5 and OpenFOAM 2D-1, the estimate of
version 1.1 and of version 1.2 (including the joint-refinement check, applied here
at every setting rather than only where the estimate is within the target) is
compared with the true error.  Reported per observable: the 10th percentile and
median of Ehat/|E_true|, and the fraction of settings where the estimate is
optimistic (below the true error), for the raw signed sum sum_k r_k and for the
final estimate of each version.  The resolution is one tenth of the tightest graded
target.
    python exp66_estimate_v12.py
"""
import os

import numpy as np

import cases as C
from errbudget import LocalBudget, V11, V121 as V12


def estimate(case, ob, a, eps, cfg, joint):
    T = case.tensor

    def look(i):
        i = tuple(i)
        if not case.usable(i) or not case.converged(i):
            return np.nan, T[i][1]
        return T[i][0][ob], T[i][1]

    lb = LocalBudget(look, case.LAD, a, hs=case.HS, fs_cap=1.25, resid_margin=cfg["resid_margin"],
                     models=case.models, kinds=case.KINDS, local_order=cfg["local_order_check"],
                     resolution=cfg["res_frac"] * eps, cancel_rule=cfg["cancel_rule"],
                     joint_res=cfg.get("joint_res", 3.0))
    lb.sweep(fine_reach=1)
    lb.per_knob_error()
    raw = float(np.nansum([r for r in lb.rk if np.isfinite(r)]))
    lb.refine_path(n_levels=0)
    if joint and np.isfinite(lb.E0):
        lb.joint_correct(look, eps, cancel_rule=cfg["cancel_rule"])
    return raw, lb.E0


def main():
    want = {"L2", "L3", "L4", "L5", "2D-1"}
    cs = [c for c in C._manufactured() if c.key in want] + [C._dfg1()]
    rows = []
    print(f"{'case':<6} {'observable':<18} {'settings':>8} | {'raw: opt.':>9} | "
          f"{'v1.1: p10':>9} {'med':>6} {'opt.':>6} | {'v1.2: p10':>9} {'med':>6} {'opt.':>6}")
    for case in cs:
        for ob, tgs in case.targets.items():
            tg = min(t for t in tgs if case.admitted(ob, t))
            eps = tg * abs(case.exact[ob])
            r_raw, r11, r12, fa11, fa12 = [], [], [], 0, 0
            for a in case.grid():
                if not (case.usable(a) and case.converged(a)):
                    continue
                Et = abs(case.tensor[a][0][ob] - case.exact[ob])
                if Et <= 0.1 * eps:          # below the resolution: no meaningful ratio
                    continue
                raw, e11 = estimate(case, ob, a, eps, V11, joint=False)
                _, e12 = estimate(case, ob, a, eps, V12, joint=True)
                r_raw.append(abs(raw) / Et); r11.append(abs(e11) / Et); r12.append(abs(e12) / Et)
                # would be accepted although outside the target: estimate within, error not
                fa11 += abs(e11) <= eps < Et
                fa12 += abs(e12) <= eps < Et
            f = lambda r: (np.percentile(r, 10), np.median(r), 100 * np.mean(np.array(r) < 1))
            a11, a12 = f(r11), f(r12)
            print(f"{case.key:<6} {ob:<18} {len(r12):>8} | {100*np.mean(np.array(r_raw) < 1):>8.1f}% | "
                  f"{a11[0]:>9.2f} {a11[1]:>6.2f} {a11[2]:>5.1f}% | {a12[0]:>9.2f} {a12[1]:>6.2f} "
                  f"{a12[2]:>5.1f}% | false acceptance at the tightest target: v1.1 {fa11}, v1.2 {fa12}", flush=True)
            rows.append(dict(case=case.key, ob=ob, n=len(r12), raw=100 * float(np.mean(np.array(r_raw) < 1)),
                             v11=a11, v12=a12, fa11=int(fa11), fa12=int(fa12), target=tg))


    import json
    json.dump(rows, open(os.path.join(C.D, "estimate_v12.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
