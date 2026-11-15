"""The three procedures of current practice (Section 6.3, S-4.2), on the stored runs.

seq        sequential convergence study (S): the parameters in the order of the case's list;
           each is refined one value at a time until the observable changes by at most eta,
           then frozen.
coord      sequential study with a joint refinement (SJ): as seq, then one run with every
           parameter refined by one step; if the observable moves by more than eta, the
           study restarts from there.
ref_check  selection against a reference run (Ref.): the cheapest setting within the target
           of the finest candidate setting.

None of the three can decline.  A run that failed or did not converge counts as failed.
"""
import os

import numpy as np


D = os.path.join(os.path.dirname(__file__), "data")


def seq(look, LAD, start, delta, order=None):
    """Sequential one-at-a-time convergence. Returns (proposal, runs)."""
    n = len(LAD)
    cur = list(start)
    runs, seen = 0, set()

    def ev(idx):
        nonlocal runs
        idx = tuple(idx)
        if idx not in seen:
            seen.add(idx)
            runs += 1
        return look(idx)[0]

    for k in (order or range(n)):
        q = ev(cur)
        while cur[k] + 1 < len(LAD[k]):
            trial = list(cur)
            trial[k] += 1
            q2 = ev(trial)
            if abs(q2 - q) <= delta:        # stopped moving: accept and freeze
                break
            cur, q = trial, q2
    return tuple(cur), runs


def coord(look, LAD, start, delta, max_pass=3):
    """SEQ plus a diagonal check, restarting if the diagonal step moves Q."""
    n = len(LAD)
    cur, runs = tuple(start), 0
    for _ in range(max_pass):
        cur, r = seq(look, LAD, cur, delta)
        runs += r
        diag = tuple(min(cur[k] + 1, len(LAD[k]) - 1) for k in range(n))
        if diag == cur:
            break
        runs += 2
        if abs(look(diag)[0] - look(cur)[0]) <= delta:
            break
        cur = diag
    return cur, runs


def ref_check(look, LAD, ALL, eps):
    """Converged-reference check: cheapest configuration within eps of the
    finest setting of every parameter. Costs the whole grid, which is what a
    reference-based study pays if it does not already own the sweep."""
    fine = tuple(len(l) - 1 for l in LAD)
    if not np.isfinite(look(fine)[0]):
        # the finest setting failed: a user would take the most refined one that ran
        fine = max(ALL, key=lambda i: (sum(i), look(i)[1]))
    qref = look(fine)[0]
    ok = [i for i in ALL if abs(look(i)[0] - qref) <= eps]
    if not ok:
        return fine, len(ALL)
    return min(ok, key=lambda i: look(i)[1]), len(ALL)


