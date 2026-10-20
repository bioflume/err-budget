"""Revision 4: the specification of Section 4 of the article, pinned to the code.

Each test states one behaviour of version 1.2.1 that the referee asked the article to
specify (look-ahead estimate, local-order check at every alignment, the exact
condition for a positive order on non-uniform candidate lists, the exponential
model, fits of nested-solver parameters, the residual gate, the joint check at its
boundaries).  If a test fails, the text of Section 4 written from these tests is
wrong; the program must not be changed to make it pass (src/errbudget is fixed by the
pre-registration of the held-out tests).

    pytest tests/test_spec_v121.py -q          # or: python tests/test_spec_v121.py
"""
import math

import numpy as np

from errbudget import LocalBudget, V121

DISC = "discretisation"
NEST = "nested-solver"


def _lb(Q, ladders, anchor, hs, kinds, eps=1.0, models=None, local=True, res_frac=0.1,
        cost=None, joint=False, cancel="none"):
    """A LocalBudget as run_protocol builds it for version 1.2.1 (resolution 0.1 eps)."""
    calls = []

    def look(i):
        i = tuple(i)
        calls.append(i)
        q = Q(i)
        c = cost(i) if cost else float(2.0 ** sum(i))
        return (q, c) if np.isfinite(q) else (np.nan, np.inf)

    lb = LocalBudget(look, ladders, anchor, hs=hs, fs_cap=1.25, resid_margin=V121["resid_margin"],
                     models=models, kinds=kinds, local_order=local, resolution=res_frac * eps,
                     cancel_rule=cancel, joint_res=V121["joint_res"])
    lb.sweep()
    lb.per_knob_error()
    lb.refine_path(n_levels=0)
    if joint and np.isfinite(lb.E0) and abs(lb.E0) <= eps:
        lb.joint_correct(look, eps, cancel_rule=cancel)
    lb.calls = calls
    return lb


# ------------------------------------------------------------------ look-ahead (Alg. 1)
def test_lookahead_is_current_estimate_plus_measured_one_step_change():
    # Phase 1 scores a move along k by |E(theta)| - |E(theta')| with
    # E(theta') = E(theta) + [Q(theta') - Q(theta)]: the current estimate (with its
    # safety factor and joint correction) plus the change measured by the sweep; no
    # parameter is re-fitted before the move is taken.
    h = [0.4, 0.2, 0.1, 0.05]
    Q = lambda i: 1.0 + 0.3 * h[i[0]] ** 2 - 0.2 * h[i[1]]
    lb = _lb(Q, [h, h], (1, 1), [h, h], [DISC, DISC], eps=1e-3)
    for k in (0, 1):
        trial = tuple(2 if j == k else 1 for j in range(2))
        assert math.isclose(lb.model_error(trial), lb.E0 + lb.e[k][2], rel_tol=0, abs_tol=1e-15)
    # the log-cost difference used in the denominator is the measured one
    assert math.isclose(lb.g[0][2] - lb.g[0][1], math.log(2.0), rel_tol=1e-12)


# --------------------------------------------------------------- local-order check
def _one_axis(qs, hs, anchor, eps=1.0, models=None, res_frac=1e-6):
    lad = [list(range(len(qs)))]
    return _lb(lambda i: qs[i[0]], lad, (anchor,), [hs], [DISC], eps=eps, models=models,
               res_frac=res_frac)


def test_positive_order_condition_is_the_exact_log_ratio_on_non_uniform_lists():
    # With h0 > h1 > h2 and changes d1, d2 of one sign, a positive order exists iff
    # d1/d2 > ln(h0/h1)/ln(h1/h2), not iff |d2| < |d1| (they agree only for uniform ratios).
    # Ratios 3/2 then 4/3 (N = 16, 24, 32): threshold 1.409.  d1/d2 = 1.2: the changes
    # decrease, yet no positive order fits them -> unbounded.
    hs = [1 / 16, 1 / 24, 1 / 32]
    d2 = 1e-2
    lb = _one_axis([0.0, 1.2 * d2, 2.2 * d2], hs, anchor=1)
    assert lb.local[0]["unbounded"]
    # Ratios 4/3 then 3/2 (N = 24, 32, 48): threshold 0.71.  d1/d2 = 0.9: the changes
    # grow, yet a positive order fits them -> bounded.
    hs = [1 / 24, 1 / 32, 1 / 48]
    lb = _one_axis([0.0, 0.9 * d2, 1.9 * d2], hs, anchor=1)
    L = lb.local[0]
    assert not L["unbounded"] and 0 < L["p"] < 30
    thr = math.log(hs[0] / hs[1]) / math.log(hs[1] / hs[2])
    assert 0.70 < thr < 0.72


def test_power_extrapolation_at_each_alignment():
    # error c h^2 on a uniform list: the check returns |c h_a^2| at the anchor, whether
    # the anchor is the middle point, the coarse end (two finer values sampled), or
    # the last point because the finer run failed.
    h = [0.4, 0.2, 0.1, 0.05]
    c = 0.5
    qs = [1.0 + c * x ** 2 for x in h]
    for a, pts in [(1, (0, 1, 2)), (0, (0, 1, 2))]:
        L = _one_axis(qs, h, anchor=a).local[0]
        assert L["pts"] == pts
        assert math.isclose(L["p"], 2.0, rel_tol=1e-6)
        assert math.isclose(L["E"], c * h[a] ** 2, rel_tol=1e-6)
    qf = qs[:3] + [np.nan]                                   # the run at h = 0.05 failed
    L = _one_axis(qf, h, anchor=2).local[0]
    assert L["pts"] == (0, 1, 2)
    assert math.isclose(L["E"], c * h[2] ** 2, rel_tol=1e-6)


def test_alternating_bound_depends_on_the_alignment():
    # changes of alternating sign that shrink: the error is bounded by the change from
    # the anchor to the last point plus 2|d2|, i.e. 3|d2| (middle), |d1 + d2| + 2|d2|
    # (coarse end) and 2|d2| (finer run failed)
    h = [0.4, 0.2, 0.1, 0.05]
    d1, d2 = 1.0e-2, -0.4e-2
    qs = [0.0, d1, d1 + d2, d1 + d2 + 1e-9]
    L = _one_axis(qs, h, anchor=1).local[0]
    assert math.isclose(L["E"], 3 * abs(d2), rel_tol=1e-9)
    L = _one_axis(qs, h, anchor=0).local[0]
    assert math.isclose(L["E"], abs(d1 + d2) + 2 * abs(d2), rel_tol=1e-9)
    L = _one_axis(qs[:3] + [np.nan], h, anchor=2).local[0]
    assert math.isclose(L["E"], 2 * abs(d2), rel_tol=1e-9)


def test_local_check_ignores_changes_below_three_resolutions():
    # The local-order check treats changes up to three times the resolution (0.3 eps
    # with the default resolution of 0.1 eps) as zero; the joint check of 1.2.1 uses
    # one resolution.  Equal changes of 0.25 eps are "converged"; of 0.35 eps,
    # unbounded (equal changes admit no positive order).
    h = [0.4, 0.2, 0.1, 0.05]
    eps = 1e-2
    for step, unbounded in ((0.25 * eps, False), (0.35 * eps, True)):
        qs = [0.0, step, 2 * step, 2 * step]
        L = _one_axis(qs, h, anchor=1, eps=eps, res_frac=0.1).local[0]
        assert L["unbounded"] == unbounded
        if not unbounded:
            assert L["E"] == 0.0


def test_exponential_axis_uses_the_mapped_variable():
    # model="exponential": x is mapped to h = exp(-(x - x_min)/span) over the
    # candidate list, so h^p = exp(-b (x - x_min)) with b = p/span; the same three-point
    # equation then gives b exactly and the error A exp(-b x_a) at the anchor.
    x = [8.0, 12.0, 16.0, 20.0]
    A, b = 2.0, 0.21
    qs = [1.0 + A * math.exp(-b * v) for v in x]
    lb = _one_axis(qs, x, anchor=1, models=["exponential"])
    L = lb.local[0]
    span = x[-1] - x[0]
    assert math.isclose(L["p"] / span, b, rel_tol=1e-6)
    assert math.isclose(L["E"], A * math.exp(-b * x[1]), rel_tol=1e-6)


# --------------------------------------------------- parameters that cannot be fitted
def test_unfittable_axis_power_versus_exponential():
    # Fewer than three successful runs (a coarse-end start without the local check):
    # a power axis is charged its one-step change times h_a/(h_a - h_b) (Eq. 4); an
    # exponential axis is charged the one-step change itself.
    h = [0.4, 0.2, 0.1]
    qs = [1.0 + 0.3 * v for v in h]
    lad = [[0, 1, 2]]
    lb = _lb(lambda i: qs[i[0]], lad, (0,), [h], [DISC], local=False)
    step = qs[1] - qs[0]
    assert math.isclose(lb.rk[0], -step * h[0] / (h[0] - h[1]), rel_tol=1e-12)
    x = [8.0, 12.0, 16.0]
    qx = [1.0 + math.exp(-0.2 * v) for v in x]
    lb = _lb(lambda i: qx[i[0]], lad, (0,), [x], [DISC], models=["exponential"], local=False)
    assert math.isclose(lb.rk[0], -(qx[1] - qx[0]), rel_tol=1e-12)


def test_nested_solver_axes_are_fitted_but_not_order_checked():
    h = [0.4, 0.2, 0.1, 0.05]
    tol = [1e-2, 1e-4, 1e-6, 1e-8]
    Q = lambda i: 1.0 + 0.1 * h[i[0]] ** 2 + 5.0 * tol[i[1]]
    lb = _lb(Q, [h, tol], (2, 2), [h, tol], [DISC, NEST], eps=1e-3)
    assert np.isfinite(lb.pk[1])            # Eq. (3) fitted along the tolerance
    assert 1 not in lb.local and 0 in lb.local


def test_residual_gate_is_the_median_misfit_of_the_fitted_axes():
    # Fs = 1.25 (1 + 5 median_k q_k) over the axes with at least three successful runs
    h = [0.4, 0.2, 0.1, 0.05]
    Q = lambda i: 1.0 + 0.1 * h[i[0]] ** 2 + 0.05 * h[i[1]] ** 1.3 + 0.01 * math.sin(7 * h[i[1]])
    lb = _lb(Q, [h, h], (3, 3), [h, h], [DISC, DISC], eps=1e-3)
    q = [x for x in lb.fit_q if np.isfinite(x)]
    assert len(q) == 2
    assert math.isclose(lb.Fs, 1.25 * (1 + 5.0 * float(np.median(q))), rel_tol=1e-12)


def test_one_step_bound_at_the_coarse_end_uses_the_second_finer_value():
    # at the coarse end the sweep samples two finer values (for the local check), and
    # the one-step bound takes the change to the finest sampled value
    h = [0.4, 0.2, 0.1, 0.05]
    qs = [1.0 + 0.5 * v ** 2 for v in h]
    lb = _one_axis(qs, h, anchor=0)
    assert lb.sampled[0] == [0, 1, 2]
    assert abs(lb.E0) >= lb.Fs * abs(lb.e[0][2]) * (1 - 1e-12)


# ----------------------------------------------------------------------- joint check
def _two_axis_problem(interaction=0.0, fail_joint=False, c1=0.01):
    """Two second-order parameters with an interaction term interaction * h0^2 h1^2.
    joint_correct runs the solver itself (the budget's own lookup), so a failed joint
    run is a failure of the solver at the jointly refined setting (3, 3)."""
    h = [0.4, 0.2, 0.1, 0.05]

    def Q(i):
        if fail_joint and tuple(i) == (3, 3):
            return np.nan
        return (1.0 + 0.02 * h[i[0]] ** 2 + c1 * h[i[1]] ** 2
                + interaction * h[i[0]] ** 2 * h[i[1]] ** 2)
    return h, Q


def test_joint_check_is_skipped_without_two_resolved_parameters():
    # only parameters whose one-step change exceeds the resolution take part; with
    # fewer than two the check is skipped and costs no run
    h = [0.4, 0.2, 0.1, 0.05]
    Q = lambda i: 1.0 + 0.02 * h[i[0]] ** 2 + 1e-9 * h[i[1]]
    lb = _lb(Q, [h, h], (2, 2), [h, h], [DISC, DISC], eps=1e-3)
    n = len(lb.calls)
    lb.joint_correct(None, 1e-3, cancel_rule="none")
    assert lb.joint["skipped"] and lb.joint["axes"] == [0] and len(lb.calls) == n


def test_joint_check_excludes_parameters_at_their_finest_value():
    h, Q = _two_axis_problem()
    Q3 = lambda i: Q(i[:2]) + 0.015 * h[i[2]] ** 2
    lb = _lb(Q3, [h, h, h], (2, 3, 2), [h, h, h], [DISC, DISC, DISC], eps=1e-3, joint=True)
    assert not lb.joint["skipped"]
    assert lb.joint["axes"] == [0, 2]
    assert lb.joint["idx"] == (3, 3, 3)


def test_joint_check_keeps_the_interaction_term_in_both_branches():
    # Ehat = Fs max(core + 2|I|, max_{k in R} |delta_k|, |delta_J|), with core = |E_fit|,
    # or sum_k |r_k| when |I| exceeds half the largest one-at-a-time change ("strong")
    # or when the signed sum cancels; the 2|I| term is kept in both branches, and the
    # estimate is only ever raised.
    eps = 1e-3
    strong = set()
    for inter in (0.5, 40.0):
        h, Q = _two_axis_problem(inter)
        lb = _lb(Q, [h, h], (2, 2), [h, h], [DISC, DISC], eps=eps)
        assert np.isfinite(lb.E0)
        before = lb.E0
        lb.joint_correct(None, eps, cancel_rule="none")
        J = lb.joint
        strong.add(bool(J["strong"]))
        rk = [r for r in lb.rk if np.isfinite(r)]
        E_fit, E_abs = sum(rk), sum(abs(r) for r in rk)
        floor = max(abs(lb.e[k][3]) for k in J["axes"])
        big = max(abs(r) for r in rk)
        cancels = abs(E_fit) < 0.8 * big
        core = (E_abs if (J["strong"] or cancels) else abs(E_fit)) + 2 * abs(J["I"])
        expect = lb.Fs * max(core, floor, abs(J["dJ"]))
        assert math.isclose(abs(J["E0_after"]), expect, rel_tol=1e-12)
        assert abs(lb.E0) == max(abs(before), abs(J["E0_after"]))     # only ever raised
    assert strong == {True, False}


def test_failed_joint_run_credits_no_cancellation():
    # if the jointly refined run fails: Ehat = Fs max(sum_k |r_k|, max_{k in R} |delta_k|)
    eps = 1e-3
    h, Q = _two_axis_problem(0.0, fail_joint=True, c1=0.016)   # both one-step changes resolved
    lb = _lb(Q, [h, h], (2, 2), [h, h], [DISC, DISC], eps=eps)
    lb.joint_correct(None, eps, cancel_rule="none")
    J = lb.joint
    assert not J["skipped"] and not np.isfinite(J["I"]) and J["dJ"] == 0.0
    rk = [r for r in lb.rk if np.isfinite(r)]
    floor = max(abs(lb.e[k][3]) for k in J["axes"])
    assert math.isclose(abs(J["E0_after"]), lb.Fs * max(sum(abs(r) for r in rk), floor), rel_tol=1e-12)


if __name__ == "__main__":
    import sys
    failed = 0
    for name, f in sorted(globals().items()):
        if name.startswith("test_") and callable(f):
            try:
                f()
                print(f"ok    {name}")
            except Exception as exc:                      # noqa: BLE001
                failed += 1
                print(f"FAIL  {name}: {type(exc).__name__}: {exc}")
    sys.exit(1 if failed else 0)
