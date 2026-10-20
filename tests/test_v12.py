"""Behaviour added in version 1.2 (revision 3)."""
import itertools

import numpy as np
import pytest

from errbudget import Knob, RunFailed, V11, V12, run_protocol, tune

KINDS2 = ["discretisation", "discretisation"]


def _replay(look, lad, eps, exact, **kw):
    """Accepted and missed starts over every start."""
    acc = miss = 0
    for a in itertools.product(*[range(len(l)) for l in lad]):
        if not np.isfinite(look(a)[0]):
            continue
        r = run_protocol(look, lad, a, eps, hs=lad, kinds=KINDS2, **kw)
        if r["certified"]:
            acc += 1
            miss += abs(look(r["proposal"])[0] - exact) > eps
    return acc, miss


def test_interacting_errors_of_opposite_sign():
    # Q = 1 + 0.3 h - 0.5 t^0.5 (1 - 2 h) (h, t ~ mesh and time step): errors of
    # opposite sign whose time-step part depends on the mesh.  Version 1.1 credits
    # the cancellation and accepts settings outside the target from every start
    # it accepts; version 1.2 credits no cancellation and accepts none that miss.
    lad = [[0.4, 0.2, 0.1, 0.05, 0.025], [0.4, 0.2, 0.1, 0.05, 0.025]]
    look = lambda i: (1.0 + 0.3 * lad[0][i[0]]
                      - 0.5 * lad[1][i[1]] ** 0.5 * (1 - 2 * lad[0][i[0]]),
                      float(2.0 ** (i[0] + i[1])))
    acc11, miss11 = _replay(look, lad, 0.03, 1.0, **V11)
    _, miss12 = _replay(look, lad, 0.03, 1.0, **V12)
    assert miss11 > 0 and miss11 == acc11
    assert miss12 == 0


def test_accelerating_convergence_at_the_coarse_end():
    # the error falls slowly over the first step and fast afterwards (a start
    # below the asymptotic range): one step from the coarse end looks converged
    x = [8.0, 12.0, 16.0, 20.0, 24.0, 28.0]
    err = [1.30, 1.15, 0.72, 0.39, 0.19, 0.11]
    lad = [x, [0.4, 0.2, 0.1]]
    look = lambda i: (1.0 + 0.01 * err[i[0]] + 1e-4 * lad[1][i[1]], float(x[i[0]] * 2 ** i[1]))
    hs = [x, lad[1]]
    for cfg, ok in ((V11, False), (V12, True)):
        r = run_protocol(look, lad, (0, 0), 0.01, hs=hs, kinds=KINDS2,
                         models=["exponential", "power"], **cfg)
        missed = r["certified"] and abs(look(r["proposal"])[0] - 1.0) > 0.01
        assert missed != ok


def test_unconverged_runs_are_failed_runs():
    tols = [1e-2, 1e-4, 1e-6, 1e-8]
    seen = []

    def solver(n, tol):
        seen.append((n, tol))
        q = 1.0 + 1.0 / n**2 + tol
        if tol == 1e-2:                  # the nested solver stops at its limit
            return {"q": q + 0.5, "converged": False}
        return {"q": q, "converged": True}

    res = tune(solver, [Knob("n", [8, 16, 32, 64]),
                        Knob("tol", tols, kind="nested-solver")],
               target=2e-3, observable="q", start={"n": 16, "tol": 1e-4})
    assert res.certified and res.settings["tol"] != 1e-2


def test_runfailed_is_recorded_not_raised():
    def solver(n):
        if n == 64:
            raise RunFailed("out of memory")
        return 1.0 + 1.0 / n**2
    res = tune(solver, [Knob("n", [8, 16, 32, 64])], target=2e-3, start={"n": 16})
    assert res.settings["n"] != 64


def test_v12_is_the_default_configuration():
    assert V12["joint_check"] and V12["local_order_check"] and V12["resid_margin"] == 5.0
    assert V12["cancel_rule"] == "none"


def test_failed_finer_run_uses_the_order_up_to_the_setting():
    # The error stalls and then falls again (an extremum that moves in jumps); the
    # run one step finer than the candidate answer fails.  Without a measurement on
    # the refining side the order is taken from the three values ending at the
    # setting: the changes grow, so the error cannot be bounded and the target is
    # declined rather than accepted on the strength of the coarse-side fit.
    err = [-77.4, -29.0, -13.1, -2.49, -2.48, -1.66, -0.61]
    x = [16, 24, 32, 48, 64, 96, 128]
    lad = [x, [1e-4, 1e-6]]
    look = lambda i: ((np.nan, 1.0) if i[0] == 6 else
                      (1.0 + 1e-4 * err[i[0]] + 1e-9 * lad[1][i[1]], float(x[i[0]] ** 2)))
    hs = [[1.0 / n for n in x], lad[1]]
    kinds = ["discretisation", "nested-solver"]
    r11 = run_protocol(look, lad, (5, 1), 1e-4, hs=hs, kinds=kinds, **V11)
    assert r11["certified"] and abs(look(r11["proposal"])[0] - 1.0) > 1e-4
    r = run_protocol(look, lad, (5, 1), 1e-4, hs=hs, kinds=kinds, **V12)
    assert not r["certified"] and r["declined"] == "unbounded_axis"


def test_v121_includes_small_resolved_changes_in_the_joint_check():
    # the frozen 1.2 skips the joint check when every one-step change is below three
    # times the resolution; 1.2.1 includes changes above the resolution itself
    from errbudget import V121
    assert V12.get("joint_res", 3.0) == 3.0 and V121["joint_res"] == 1.0
    assert {k: v for k, v in V121.items() if k != "joint_res"} == V12
