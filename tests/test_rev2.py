"""Behaviour added in version 1.1 (revision 2)."""
import itertools

import numpy as np
import pytest

from knobs import Knob, measure_floor, run_protocol, tune


def _sep(ladders, p=(2.0, 1.0)):
    """Separable synthetic solver: Q = 1 + sum_k h_k^p_k, cost grows as refined."""
    def look(idx):
        q = 1.0 + sum(ladders[k][j] ** p[k] for k, j in enumerate(idx))
        c = float(np.prod([2.0 ** j for j in idx]))
        return q, c
    return look


def test_noise_floor_declines_without_running():
    lad = [[0.4, 0.2, 0.1, 0.05], [0.4, 0.2, 0.1, 0.05]]
    calls = []
    look = _sep(lad)
    r = run_protocol(lambda i: (calls.append(i), look(i))[1], lad, (1, 1), 1e-3,
                     hs=lad, noise_floor=1e-3)
    assert not r["certified"] and r["declined"] == "below_floor" and not calls


def test_measure_floor_is_the_spread():
    assert measure_floor([1.0, 1.5, 0.75]) == pytest.approx(0.75)
    with pytest.raises(ValueError):
        measure_floor([1.0])


def test_failed_finer_run_is_not_treated_as_converged():
    lad = [[0.4, 0.2, 0.1, 0.05], [0.4, 0.2, 0.1, 0.05]]
    good = _sep(lad)
    # the finest value of the second parameter always fails (e.g. out of memory)
    look = lambda i: (np.nan, np.inf) if i[1] == 3 else good(i)
    exact = 1.0
    for a in itertools.product(range(4), range(3)):
        r = run_protocol(look, lad, a, 0.03, hs=lad)
        if r["certified"]:
            assert abs(look(r["proposal"])[0] - exact) <= 0.03


def test_exponential_knob_requires_increasing_values():
    Knob("ecut", [10, 20, 30], model="exponential")
    with pytest.raises(ValueError):
        Knob("ecut", [30, 20, 10], model="exponential")
    with pytest.raises(ValueError):
        Knob("ecut", [10, 20, 30], model="gaussian")


def test_exponential_model_recovers_exponential_convergence():
    ec = [8.0, 12.0, 16.0, 20.0, 24.0, 28.0]
    def solver(ecut, tol):
        return 1.0 + 0.5 * np.exp(-0.4 * ecut) + tol, ecut ** 1.5
    res = tune(solver, [Knob("ecut", ec, model="exponential"),
                        Knob("tol", [1e-2, 1e-4, 1e-6, 1e-8], kind="nested-solver")],
               target=1e-4, start={"ecut": 16.0, "tol": 1e-4}, screen=False)
    assert res.certified
    q, _ = solver(**res.settings)
    assert abs(q - 1.0) <= 1e-4


def test_tune_floor_probe_declines():
    calls = {"n": 0}
    def solver(n):
        calls["n"] += 1
        return 1.0 + 1.0 / n**2 + 1e-3 * (calls["n"] % 2), 1.0 * n
    res = tune(solver, [Knob("n", [4, 8, 16, 32])], target=1e-4,
               floor_probe=[{"n": 64}, {"n": 64}, {"n": 128}], screen=False)
    assert not res.certified
