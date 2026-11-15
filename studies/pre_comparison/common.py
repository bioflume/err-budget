"""Shared set-up: the OpenFOAM 2D-1 study, its graded test cases, errbudget's replay and the
PRE selector.  The study, the program and the grading are read from the article's
repository (read-only); nothing there is modified.

    ERRBUDGET_REPO   path of the article's repository (default: the repository holding this folder)
"""
from __future__ import annotations

import os
import sys

import numpy as np

REPO = os.environ.get("ERRBUDGET_REPO",
                      os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.dont_write_bytecode = True          # write no __pycache__ into the read-only repository
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "paper"))

import cases as C                                   # noqa: E402  (the article's case registry)
import pooled as P                                  # noqa: E402  (its configurations)
from errbudget import run_protocol                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
CASE = C._dfg1()                                    # cost = mesh cells x outer iterations
CFG = P.CONFIGS["v1.2.1"]                           # the released program, as in Table S-3.1
NC_DEFAULT = 0                                      # nCellsInCoarsestLevel = 10, OpenFOAM's default
TESTS = [("cd", 1e-2), ("cd", 3e-3), ("cd", 1e-3), ("dp", 3e-2), ("dp", 1e-2),
         ("cl", 1e-1), ("cl", 3e-2)]
ZS = (1.0, 2.0, 3.0)


def usable(i):
    return CASE.usable(i)


def good(i):
    """A run that completed and met its own convergence criterion."""
    return CASE.usable(i) and CASE.converged(i)


def cost(i):
    return float(CASE.tensor[tuple(i)][1])


def value(i, ob):
    return float(CASE.tensor[tuple(i)][0][ob])


def x_of(i):
    """PRE's inputs: errbudget's refinement scalars of the mesh and the two tolerances,
    h = (1 / n_arc, residualControl, GAMG tolerance); the cost-only parameter is dropped."""
    return (CASE.HS[0][i[0]], CASE.HS[1][i[1]], CASE.HS[2][i[2]])


def eps_abs(ob, tg):
    return tg * abs(CASE.exact[ob])


def optimum(ob, tg, grid=None):
    """Cheapest converged setting meeting the target (the denominator of the cost ratio)."""
    eps = eps_abs(ob, tg)
    G = grid if grid is not None else CASE.grid()
    feas = [i for i in G if good(i) and abs(value(i, ob) - CASE.exact[ob]) <= eps]
    return min(feas, key=cost)


def starts():
    """errbudget's graded starts: every converged setting (608)."""
    return [i for i in CASE.grid() if good(i)]


def errbudget_search(ob, tg, a, ladders=None, tensor_index=None):
    """One search of the released program from start a, exactly as exp72_baselines_pooled.cell
    replays it: unconverged and failed runs are paid for and returned as failed.  Returns the
    program's result and the runs it used, in the order it asked for them.

    ladders / tensor_index: optional restriction (e.g. the 3-parameter problem at the
    default nCellsInCoarsestLevel); tensor_index maps the restricted index to the full one."""
    eps = eps_abs(ob, tg)
    T = CASE.tensor
    LAD = CASE.LAD if ladders is None else ladders
    full = (lambda i: tuple(i)) if tensor_index is None else tensor_index
    seen = {}

    def look(i):
        i = tuple(i)
        if i not in seen:
            j = full(i)
            if not good(j):
                seen[i] = (np.nan, (T[j][1] if np.isfinite(T[j][1]) else np.inf))
            else:
                seen[i] = (T[j][0][ob], T[j][1])
        return seen[i]

    kw = dict(CFG["kw"])
    if ladders is None:
        hs, kinds = CASE.HS, list(CASE.KINDS)
    else:
        hs, kinds = CASE.HS[:len(LAD)], list(CASE.KINDS[:len(LAD)])
    r = run_protocol(look, LAD, a, eps, hs=hs, kinds=kinds, models=CASE.models, **kw)
    used = [full(i) for i in seen]
    return r, used


def pre_dataset(runs, ob):
    """PRE's training data from a set of runs: converged runs only, keyed by the PRE input x;
    runs that differ only in the cost-only parameter are averaged (their observables agree to
    about 1e-7 relative).  Returns X (n, 3), Y (n,), and the runs behind each row."""
    rows = {}
    for i in runs:
        if not good(i):
            continue
        rows.setdefault(x_of(i), []).append(i)
    keys = sorted(rows)
    X = np.array(keys, dtype=float)
    Y = np.array([np.mean([value(i, ob) for i in rows[k]]) for k in keys])
    return X, Y, [rows[k] for k in keys]


def select(runs, ob, tg, m, s, z):
    """The protocol's selector: the cheapest converged run among those PRE has seen with
    |Q(theta) - m| + z s <= eps; None (decline) if no run qualifies or the fit failed."""
    if not (np.isfinite(m) and np.isfinite(s)):
        return None
    eps = eps_abs(ob, tg)
    ok = [i for i in runs if good(i) and abs(value(i, ob) - m) + z * s <= eps]
    return min(ok, key=cost) if ok else None


def grade(q, ob, tg):
    """True relative error of a returned setting over the target (> 1: missed)."""
    return abs(value(q, ob) - CASE.exact[ob]) / eps_abs(ob, tg)
