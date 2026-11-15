"""Convergence thresholds chosen once (S-4.1, Table S-4.1): the failure rate defined in
exp47_transfer.py, without and with the two adjustments.  With them, deviations and
remaining errors below a tenth of the tightest graded target count as zero, and runs whose
iterative solver did not converge count as failed.

    python exp64_transfer_floor.py      # writes data/transfer_floor.json
"""
import itertools
import os

import numpy as np

import cases as C

PHIS = [0.01, 0.05, 0.1, 0.25, 0.5]
NAMES = {"L2": ["N", "time steps", "Newton tol", "GMRES tol"],
         "L3": ["N", "time steps", "Newton tol", "GMRES tol", "ILU drop"],
         "L4": ["N", "split steps", "react substeps", "GMRES tol", "ILU drop"],
         "L5": ["N", "Newton tol", "GMRES tol", "ILU drop"],
         "2D-1": ["n_arc", "residualControl", "GAMG tol", "coarsest level"],
         "2D-2": ["n_arc", "dt", "residualControl", "GAMG tol"],
         "QE": ["ecutwfc", "k-mesh", "conv_thr", "mixing_beta"],
         "ABINIT": ["ecut", "ngkpt", "toldfe", "nline"]}


def failure_rate(case, k, ob, phi, floor, fail_unconverged=True):
    LAD, T = case.LAD, case.tensor
    others = [j for j in range(len(LAD)) if j != k]
    combos = list(itertools.product(*[range(len(LAD[j])) for j in others]))
    fin = tuple(len(LAD[j]) - 1 for j in others)
    ex = case.exact[ob]

    def val(idx):
        if not case.usable(idx) or (fail_unconverged and not case.converged(idx)):
            return np.nan
        return T[idx][0][ob]

    def thr(c):
        idx = [0] * len(LAD)
        for kk, v in zip(others, c):
            idx[kk] = v
        vals = []
        for j in range(len(LAD[k])):
            idx[k] = j
            vals.append(val(tuple(idx)))
        if not all(np.isfinite(v) for v in vals):
            return None
        fine = vals[-1]
        tot = abs(fine - ex)
        tol = max(phi * tot, floor)
        for j, v in enumerate(vals):
            if abs(v - fine) <= tol:
                return j
        return len(LAD[k]) - 1

    ref = thr(fin)
    if ref is None:                    # fall back to the most refined complete sweep
        for c in sorted(combos, key=lambda c: (-sum(c), [-x for x in c])):
            ref = thr(c)
            if ref is not None:
                break
    if ref is None:
        return 0, 0
    n = fail = 0
    for c in combos:
        t = thr(c)
        if t is None:
            continue
        n += 1
        fail += t > ref
    return fail, n


def main():
    want = {"L2", "L3", "L4", "L5", "2D-1", "2D-2", "QE", "ABINIT"}
    rows = []
    for case in C.all_cases():
        if case.key not in want:
            continue
        obs = [o for o in case.targets if any(case.admitted(o, t) for t in case.targets[o])]
        print(f"\n{case.label}")
        print(f"  {'parameter':<16} {'class':<15} {'floor':<6} " + " ".join(f"phi={p:<5g}" for p in PHIS))
        for k in range(len(case.LAD)):
            for use_floor in (False, True):
                cells = []
                for phi in PHIS:
                    F = N = 0
                    for ob in obs:
                        tmin = min(t for t in case.targets[ob] if case.admitted(ob, t))
                        fl = 0.1 * tmin * abs(case.exact[ob]) if use_floor else 0.0
                        f, n = failure_rate(case, k, ob, phi, fl, fail_unconverged=use_floor)
                        F += f; N += n
                    cells.append(100.0 * F / max(N, 1))
                rows.append(dict(case=case.key, param=NAMES[case.key][k], kind=case.KINDS[k],
                                 floor=use_floor, rates=cells))
                print(f"  {NAMES[case.key][k]:<16} {case.KINDS[k]:<15} {'yes' if use_floor else 'no':<6} "
                      + " ".join(f"{x:>8.0f}%" for x in cells), flush=True)
    import json
    json.dump(rows, open(os.path.join(C.D, "transfer_floor.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
