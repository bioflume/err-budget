"""What the method is actually competing with.

Three baselines, each a written-down version of what practitioners do, graded on
the same tensors and from the same starting configurations as the method.

  SEQ   Sequential one-at-a-time convergence. Take the parameters in the order
        given. For each, refine one rung at a time until the observable stops
        moving by more than a threshold delta, then freeze it and move on. This
        is what every convergence workchain in the workflow frameworks does, and
        what a careful user does by hand.

  COORD Coordinate search with a diagonal check, after Grossmann et al.: run SEQ,
        then verify by advancing every parameter one rung together; if the
        observable moves by more than delta, restart from that point. Costs more
        and catches the case where the parameters are not separable.

  REF   Converged-reference check, after Khouzami et al.: compute a reference at
        the finest setting of every parameter, then accept the cheapest
        configuration on the grid whose observable is within the target of that
        reference. This is the strongest baseline -- it uses the whole grid ---
        and it is also the one whose assumption a discretised PDE violates,
        because the reference is not the truth.

All three need a threshold. A user has no principled way to pick one, so each is
run at delta = eps, eps/3 and eps/10 relative to the target and the best and
worst outcomes are both reported. That is the fairest treatment: the method gets
no threshold to tune, and the baselines get three.
"""
import itertools
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from knobs import run_protocol
from knobs.problems import Burgers1D, Burgers2D, ReactionDiffusion2D

D = os.path.join(os.path.dirname(__file__), "data")


# ------------------------------------------------------------------ baselines
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



