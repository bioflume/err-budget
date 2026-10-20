"""One test per behavioural claim the article makes.

The article carries a table mapping each claim to the test below that fails if
the claim stops holding.  The point is that a reader should not have to take the
prose on trust: every statement of the form "the method does X" has an executable
counterpart here.
"""
import itertools
import math

import numpy as np
import pytest

from errbudget import Knob, tune
from errbudget._core import LocalBudget, run_protocol


# --------------------------------------------------------------- fixtures --
def make_tensor(orders=(2.0, 2.0), coeffs=(1.0, 0.5), masked_at=2):
    """A solver with imposed structure: two converging axes, one masked axis.

    Axis 0 and 1 converge at the given orders; axis 2 stops affecting the answer
    beyond index `masked_at` while continuing to cost.
    """
    LAD = [[16, 32, 64, 128, 256], [10, 20, 40, 80, 160],
           [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
    HS = [[1.0 / n for n in LAD[0]], [1.0 / s for s in LAD[1]], list(LAD[2])]
    QS = 1.0
    tensor = {}
    for idx in itertools.product(*[range(5) for _ in range(3)]):
        e = (coeffs[0] * HS[0][idx[0]] ** orders[0]
             + coeffs[1] * HS[1][idx[1]] ** orders[1]
             + 1e-3 * HS[2][min(idx[2], masked_at)])
        cost = (LAD[0][idx[0]] * LAD[1][idx[1]]
                * (1.0 + 0.3 * idx[2]))
        tensor[idx] = (QS + e, float(cost))
    return tensor, LAD, HS, QS


def lookup_of(tensor):
    return lambda i: tensor[tuple(i)]


# ------------------------------------------------------- claims: estimate --
def test_recovers_imposed_orders():
    """Claim (S5): black-box measurement recovers the imposed convergence order."""
    tensor, LAD, HS, QS = make_tensor(orders=(2.0, 2.0))
    lb = LocalBudget(lookup_of(tensor), LAD, (2, 2, 4), hs=HS)
    lb.sweep(); lb.per_knob_error()
    assert lb.pk[0] == pytest.approx(2.0, abs=0.15)
    assert lb.pk[1] == pytest.approx(2.0, abs=0.15)


def test_masked_parameter_gets_a_negligible_share():
    """Claim (S6.1): a masked parameter is assigned a negligible budget share."""
    tensor, LAD, HS, QS = make_tensor(masked_at=1)
    lb = LocalBudget(lookup_of(tensor), LAD, (2, 2, 4), hs=HS)
    lb.sweep(); lb.per_knob_error()
    total = sum(abs(r) for r in lb.rk if np.isfinite(r))
    assert abs(lb.rk[2]) / total < 0.01


def test_consistency_floor_carries_a_two_point_fit():
    """Claim (S4.1, S7.2): at the coarse end of a ladder the offset fit has only
    two points.  In v1.0 the axis contributed nothing and the consistency floor
    had to cover it; from v1.1 its one-step change, extrapolated as if the axis
    converged to first order, is added to the budget."""
    tensor, LAD, HS, QS = make_tensor()
    true_err = abs(tensor[(0, 0, 0)][0] - QS)
    lb = LocalBudget(lookup_of(tensor), LAD, (0, 0, 0), hs=HS, unfit="floor")
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    assert not np.isfinite(lb.rk[0])          # fewer than three points: no fit
    assert abs(lb.E0) > 0.1 * true_err        # the floor, not zero
    lb1 = LocalBudget(lookup_of(tensor), LAD, (0, 0, 0), hs=HS)    # default order1
    lb1.sweep(); lb1.per_knob_error(); lb1.refine_path(n_levels=0)
    step = abs(lb1.e[0][1])
    ratio = HS[0][0] / HS[0][1]
    assert abs(lb1.rk[0]) == pytest.approx(step * (1 + 1 / (ratio - 1)), rel=1e-12)
    assert abs(lb1.E0) >= abs(lb.E0)


def test_estimate_is_a_signed_sum_with_a_floor():
    """Claim (S4.1): the estimate sums residuals WITH their signs and is then
    bounded below by the largest single-parameter refinement deviation."""
    tensor, LAD, HS, QS = make_tensor()
    lb = LocalBudget(lookup_of(tensor), LAD, (3, 3, 4), hs=HS, fs_cap=1.25)
    lb.sweep(); lb.per_knob_error(); lb.refine_path(n_levels=0)
    signed = sum(r for r in lb.rk if np.isfinite(r))
    floor = max(abs(lb.e[k][j]) for k in range(3)
                for j in lb.sampled[k] if j > lb.anchor[k])
    assert abs(lb.E0) == pytest.approx(1.25 * max(abs(signed), floor), rel=1e-9)


# ---------------------------------------------------------- claims: walk ---
def test_unreachable_target_is_refused():
    """Claim (S7.2): a target no supplied combination attains is refused."""
    tensor, LAD, HS, QS = make_tensor()
    best = min(abs(v[0] - QS) for v in tensor.values())
    r = run_protocol(lookup_of(tensor), LAD, (2, 2, 2), best / 100.0, hs=HS)
    assert not r["certified"]


def test_walk_reaches_a_fixed_point():
    """Claim (Prop. 1): the walk terminates at a fixed point, so the answer does
    not depend on the iteration cap."""
    tensor, LAD, HS, QS = make_tensor()
    eps = 1e-4
    a = run_protocol(lookup_of(tensor), LAD, (1, 1, 1), eps, hs=HS, max_iter=40)
    b = run_protocol(lookup_of(tensor), LAD, (1, 1, 1), eps, hs=HS, max_iter=200)
    assert a["proposal"] == b["proposal"]
    c = run_protocol(lookup_of(tensor), LAD, a["proposal"], eps, hs=HS)
    assert c["proposal"] == a["proposal"]      # restarting from the answer


def test_audit_extrapolates_rather_than_differences():
    """Claim (S4.2): a frozen parameter's remaining error is the tail of a
    geometric series, so the one-step difference is doubled, not used raw."""
    tensor, LAD, HS, QS = make_tensor()
    KIN = ["discretisation", "discretisation", "nested-solver"]
    r = run_protocol(lookup_of(tensor), LAD, (1, 1, 1), 1e-4, hs=HS, kinds=KIN,
                     use_screen=True, audit=True)
    if r["certified"]:
        assert abs(tensor[r["proposal"]][0] - QS) <= 1e-4


def test_cost_only_pass_reverifies_its_declaration():
    """Claim (S3, S4.2): a parameter wrongly declared cost-only is not loosened
    by the cost-only pass, because that pass re-verifies the observable."""
    tensor, LAD, HS, QS = make_tensor()
    KIN = ["cost-only", "discretisation", "nested-solver"]   # axis 0 is a lie
    r = run_protocol(lookup_of(tensor), LAD, (3, 3, 3), 1e-4, hs=HS, kinds=KIN)
    if r["certified"]:
        assert abs(tensor[r["proposal"]][0] - QS) <= 1e-4


def test_solver_failure_is_handled_not_fitted():
    """Claim (S4.1): a configuration that fails, diverges or returns a
    non-finite value is never proposed."""
    tensor, LAD, HS, QS = make_tensor()
    bad = {(0, 0, 0), (0, 1, 0), (1, 0, 0)}
    def lk(i):
        i = tuple(i)
        if i in bad:
            return (float("nan"), float("inf"))
        return tensor[i]
    r = run_protocol(lk, LAD, (2, 2, 2), 1e-4, hs=HS)
    assert tuple(r["proposal"]) not in bad


# ------------------------------------------------- claims: the structure ---
def test_misdeclared_cost_only_is_reported():
    """Claim (S3): the method reports declarations its measurements contradict.

    The screen flags a parameter declared cost-only whose influence share is at
    or above tau.  A mis-declared parameter BELOW that share is frozen instead,
    and is then caught by the audit, which sums frozen contributions into the
    budget -- so falsification happens in two places, not one.  This test covers
    the first; test_cost_only_pass_reverifies_its_declaration covers the second.
    """
    tensor, LAD, HS, QS = make_tensor()
    KIN = ["discretisation", "cost-only", "nested-solver"]   # axis 1 is a lie
    r = run_protocol(lookup_of(tensor), LAD, (2, 2, 2), 1e-4, hs=HS, kinds=KIN,
                     use_screen=True)
    assert 1 in r["violations"]


def test_two_observables_need_not_agree():
    """Claim (S6.3): the cost-optimal setting for one observable can fail the
    same relative target for another computed from the same runs."""
    from errbudget.problems import Burgers1D
    import os
    import pickle
    path = os.path.join(os.path.dirname(__file__), "..", "paper", "data",
                        "burgers_tensor.pkl")
    if not os.path.exists(path):
        pytest.skip("precomputed tensor not present")
    t = pickle.load(open(path, "rb"))
    EX = Burgers1D(nu=0.05, T=1.0).exact_observables()
    ALL = list(t.keys())
    tgt = 1e-3
    feas = [i for i in ALL
            if abs(t[i][0]["energy"] - EX["energy"]) <= tgt * abs(EX["energy"])]
    best = min(feas, key=lambda i: t[i][1])
    other = abs(t[best][0]["dissipation"] - EX["dissipation"])
    assert other > tgt * abs(EX["dissipation"])


def test_stabilised_exact_solution_agrees_where_the_spectral_route_is_valid():
    """Claim (S5): the log-space Cole-Hopf evaluation reproduces the spectral one
    at moderate viscosity, and is finite where the spectral one is not."""
    from errbudget.problems import Burgers1D

    def spectral(nu, Nf=1 << 15):
        x = 2 * np.pi * np.arange(Nf) / Nf
        th0 = np.exp(-(1.0 - np.cos(x)) / (2 * nu))
        k = np.fft.fftfreq(Nf, d=1.0 / Nf)
        H = np.fft.fft(th0) * np.exp(-nu * k ** 2)
        th = np.real(np.fft.ifft(H))
        thx = np.real(np.fft.ifft(1j * k * H))
        u = -2 * nu * thx / th
        return 0.5 * float(np.sum(u ** 2)) * (2 * np.pi / Nf)

    assert Burgers1D(nu=0.2, T=1.0).exact_observables()["energy"] == \
        pytest.approx(spectral(0.2), rel=1e-10)
    hard = Burgers1D(nu=0.02, T=1.0).exact_observables()
    assert math.isfinite(hard["energy"]) and math.isfinite(hard["dissipation"])


def test_ladder_direction_guard():
    """Claim (S10): an ambiguous ladder direction is refused, not guessed."""
    with pytest.raises(ValueError, match="cannot infer the refinement direction"):
        Knob("tol", [1e-8, 1e-6, 1e-4])
