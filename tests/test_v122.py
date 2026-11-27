"""Engine 1.2.2 and errbudget 2.1: one test per correction (CHANGELOG.md), and tests for
behaviours that planted faults showed no other test checked.

Written by a coding agent; reviewed and accepted by the author on 27 November 2026, and
listed in PROTECTED since (AGENTS.md).

    pytest tests/test_v122.py -q
"""
import json
import math
import sys

import numpy as np
import pytest

from errbudget import (LocalBudget, Parameter, RunFailed, V121, V122, full_configuration,
                       run_protocol, tune)
from errbudget import api as api_module
from errbudget.jobfile import JobSpec, load_job
from errbudget.runner import CommandSolver, JobError, extract_number, prepare_run_dir

DISC, NEST = "discretisation", "nested-solver"


def _table(Q, ladders, cost=None, fail=()):
    """lookup over a full grid: Q(idx) -> value; failed settings give (nan, inf)."""
    def look(i):
        i = tuple(i)
        if i in fail:
            return np.nan, np.inf
        return Q(i), (cost(i) if cost else float(2.0 ** sum(i)))
    return look


def _lb(Q, ladders, anchor, hs, kinds, eps, cfg, cost=None):
    """A LocalBudget built and run as run_protocol's analyse() builds it."""
    look = _table(Q, ladders, cost)
    lb = LocalBudget(look, ladders, anchor, hs=hs, fs_cap=1.25, resid_margin=cfg["resid_margin"],
                     kinds=kinds, local_order=cfg["local_order_check"],
                     resolution=cfg["res_frac"] * eps, cancel_rule=cfg["cancel_rule"],
                     joint_res=cfg.get("joint_res", 3.0),
                     strict_unfit=cfg.get("strict_unfit", False),
                     local_at_finest=cfg.get("local_at_finest", False),
                     cancel_strict=cfg.get("cancel_strict", False))
    lb.sweep()
    lb.per_knob_error()
    lb.refine_path(n_levels=0)
    return lb


# ------------------------------------------------------------- the three corrections
def test_short_axis_with_a_failed_neighbour_is_not_error_free():
    # three values, the middle run fails, start at the finest: 1.2.1 certified it with
    # 98 times the target error; 1.2.2 charges the coarse change extrapolated to first
    # order, which bounds the true error
    N = [8, 16, 32]
    look = _table(lambda i: 1 + 1 / N[i[0]] ** 2, [N], fail={(1,)})
    kw = dict(hs=[[1 / n for n in N]], kinds=[DISC], models=["power"])
    old = run_protocol(look, [N], (2,), 1e-5, **kw, **V121)
    new = run_protocol(look, [N], (2,), 1e-5, **kw, **V122)
    assert old["certified"] and math.isnan(old["budget"].rk[0])
    assert not new["certified"]
    assert new["budget"].rk[0] >= 1 / 32 ** 2          # at least the true error


def test_swept_axis_with_no_finite_change_is_unmeasured():
    N = [8, 16, 32]
    look = _table(lambda i: 1 + 1 / N[i[0]] ** 2, [N], fail={(0,), (1,)})
    r = run_protocol(look, [N], (2,), 1e-5, hs=[[1 / n for n in N]], kinds=[DISC],
                     models=["power"], **V122)
    assert not r["certified"] and r["declined"] == "unmeasured_axis"


def test_local_order_check_runs_at_the_finest_value():
    h = [0.4, 0.2, 0.1, 0.05]
    Q = lambda i: 1 + h[i[0]] ** 2 - 4 * h[i[0]] ** 3
    for cfg, expect in ((V121, False), (V122, True)):
        lb = _lb(Q, [h], (3,), [h], [DISC], 1e-6, cfg)
        assert (0 in lb.local) is expect
    lb = _lb(Q, [h], (3,), [h], [DISC], 1e-6, V122)
    assert lb.local[0]["pts"] == (1, 2, 3)


def test_no_credit_for_any_cancellation():
    # three axes whose errors at the anchor are +1, +1, -1 (in units of h_a^2): 1.2.1
    # credited the cancellation (net 1 is above 0.8 times the largest term), 1.2.2 sums
    # the magnitudes
    h = [0.4, 0.2, 0.1, 0.05]
    coef = (1.0, 1.0, -1.0)
    Q = lambda i: 1 + sum(c * h[j] ** 2 for c, j in zip(coef, i))
    kinds = [DISC] * 3
    old = _lb(Q, [h] * 3, (2, 2, 2), [h] * 3, kinds, 1.0, V121)
    new = _lb(Q, [h] * 3, (2, 2, 2), [h] * 3, kinds, 1.0, V122)
    total = sum(abs(r) for r in new.rk)
    assert not old.cancelled and new.cancelled
    assert abs(new.E0) >= total * 1.25 * (1 - 1e-9)
    assert abs(old.E0) < total


def test_v121_configuration_is_unchanged():
    assert {k: v for k, v in V122.items()
            if k not in ("strict_unfit", "local_at_finest", "cancel_strict")} == V121
    full = full_configuration("v1.2.1")
    assert not full["strict_unfit"] and not full["local_at_finest"] and not full["cancel_strict"]


# ------------------------------------------------------------------- configurations
def test_named_v12_runs_with_its_own_joint_resolution(monkeypatch):
    seen = {}
    real = api_module.run_protocol

    def spy(*a, **k):
        seen.update(k)
        return real(*a, **k)

    monkeypatch.setattr(api_module, "run_protocol", spy)
    tune(lambda n: 1 + 1 / n ** 2, [Parameter("n", [8, 16, 32, 64, 128])], 1e-3,
         configuration="v1.2")
    assert seen["joint_res"] == 3.0 and seen["strict_unfit"] is False
    tune(lambda n: 1 + 1 / n ** 2, [Parameter("n", [8, 16, 32, 64, 128])], 1e-3)
    assert seen["joint_res"] == 1.0 and seen["strict_unfit"] is True


def test_job_file_configuration_v11_is_v11():
    d = {"command": "x {{n}}", "parameters": [{"name": "n", "values": [1, 2]}],
         "observable": "v", "target": 1, "options": {"configuration": "v1.1"}}
    kw = JobSpec.from_dict(d).tune_kwargs()
    assert kw["joint_res"] == 3.0 and kw["strict_unfit"] is False and kw["joint_check"] is False


# ------------------------------------------------------------------- tune
def test_relative_target_holds_for_the_converged_value():
    # a coarse start overstates |Q| by 78%; the relative target must still hold for Q*
    calls = []

    def solver(N):
        calls.append(N)
        return 1 + 50 / N ** 2

    N = [8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    r = tune(solver, [Parameter("N", N)], 1e-3, relative=True, start={"N": 8})
    if r.certified:
        q = 1 + 50 / r.settings["N"] ** 2
        assert abs(q - 1.0) <= 1e-3 * 1.0
    assert len(calls) == len(set(calls))               # every setting run once


def test_relative_target_with_a_failed_start_says_so():
    def solver(n):
        if n == 32:
            raise RunFailed("crash")
        return 1 + 1 / n ** 2

    with pytest.raises(ValueError, match="relative target"):
        tune(solver, [Parameter("n", [8, 16, 32, 64])], 1e-3, relative=True)


def test_declined_result_says_why():
    def solver(n):
        if n == 32:
            raise RunFailed("crash")
        return 1 + 1 / n ** 2

    r = tune(solver, [Parameter("n", [8, 16, 32, 64])], 1e-3)
    assert not r.certified and "start setting" in r.declined
    assert "Why:" in str(r) and "Why:" in r.markdown()


def test_floor_probe_survives_a_failed_probe_and_refuses_too_few():
    def solver(n, fail=()):
        if n in fail:
            raise RunFailed("out of memory")
        return 1 + 1 / n ** 2 + 1e-9 * (n % 3)

    p = [Parameter("n", [8, 16, 32, 64])]
    r = tune(lambda n: solver(n, fail=(1024,)), p, 1e-3,
             floor_probe=[{"n": 256}, {"n": 512}, {"n": 1024}])
    assert r.certified
    with pytest.raises(ValueError, match="floor_probe"):
        tune(lambda n: solver(n, fail=(512, 1024)), p, 1e-3,
             floor_probe=[{"n": 512}, {"n": 1024}])
    with pytest.raises(ValueError, match="give every parameter"):
        tune(solver, p, 1e-3, floor_probe=[{"m": 1}, {"m": 2}])


def test_failed_run_cost_is_the_reported_one():
    def solver(n):
        if n == 64:
            raise RunFailed("crash", cost=1000.0)
        return 1 + 1 / n ** 2, float(n)

    r = tune(solver, [Parameter("n", [8, 16, 32, 64])], 1e-8)
    assert r.calibration_cost >= 1000.0


def test_unconverged_runs_are_failed_by_default():
    def solver(n):
        return {"u": 1 + 1 / n ** 2 + (1.0 if n == 64 else 0.0), "converged": n != 64}

    p = [Parameter("n", [8, 16, 32, 64, 128, 256])]
    r = tune(solver, p, 1e-4, observable="u", start={"n": 32})
    assert r.settings["n"] != 64
    r2 = tune(solver, p, 1e-4, observable="u", start={"n": 32}, accept_unconverged=True)
    assert r2.settings != r.settings or not r2.certified or r2.error_estimate != r.error_estimate


def test_markdown_shows_values_in_full_and_terminal_shows_frozen():
    from errbudget import TuneResult
    r = TuneResult(settings={"a": 1.2345678}, certified=True, error_estimate=1e-4,
                   target=1e-3, observable="q", limiting=None, budget_shares={}, wasted=[],
                   frozen=["b"], violations=[], runs_used=3, calibration_cost=3.0,
                   recommended_cost=1.0, start={"a": 1.2345679})
    assert "1.2345678" in r.markdown() and "1.2345679" in r.markdown()
    assert "frozen by the preliminary sweep" in str(r)


# -------------------------------------------------------- behaviours that planted faults
# changed without any test noticing (TRUST.md, "Planted faults")
def test_fitted_remaining_error_has_the_right_sign_and_size():
    # e(h) = c h^2 with c > 0: the anchor carries +c h_a^2 (a flipped fit basis gives -)
    h = [0.4, 0.2, 0.1, 0.05, 0.025]
    for c in (2.0, -3.0):
        lb = _lb(lambda i: 1 + c * h[i[0]] ** 2, [h], (3,), [h], [DISC], 1.0, V122)
        # the order is fitted on a grid (p = 2.02 here), so the size agrees to a few per cent
        assert np.sign(lb.rk[0]) == np.sign(c)
        assert lb.rk[0] == pytest.approx(c * h[3] ** 2, rel=0.1)


def test_fit_quality_is_the_relative_rms_residual():
    # fit_q is sqrt(residual / sum e^2), which the residual gate scales the margin by
    h = [0.4, 0.2, 0.1, 0.05, 0.025]
    Q = lambda i: 1 + h[i[0]] ** 2 + 0.01 * (-1) ** i[0]
    lb = _lb(Q, [h], (4,), [h], [DISC], 1.0, V121)
    p, ha = lb.pk[0], h[4]
    js = [j for j in lb.sampled[0] if np.isfinite(lb.e[0][j])]
    e = lb.e[0][js]
    c = lb.rk[0] / ha ** p
    fit = c * (np.asarray(h)[js] ** p - ha ** p)
    assert lb.fit_q[0] == pytest.approx(np.sqrt(np.sum((e - fit) ** 2) / np.sum(e ** 2)), rel=1e-9)


def test_protocol_cost_counts_each_run_once():
    h = [0.4, 0.2, 0.1, 0.05, 0.025]
    costs = []

    def look(i):
        c = float(3 ** i[0] * 2 ** i[1])
        costs.append((tuple(i), c))
        return 1 + h[i[0]] ** 2 + 0.1 * h[i[1]], c

    r = run_protocol(look, [h, h], (2, 2), 1e-4, hs=[h, h], kinds=[DISC, DISC],
                     models=["power"] * 2, **V122)
    assert r["protocol_cost"] == pytest.approx(sum(c for _, c in dict(costs).items()))
    assert r["runs"] == len(dict(costs))


def test_phase_one_refines_the_parameter_that_holds_the_error():
    h = [0.4, 0.2, 0.1, 0.05, 0.025, 0.0125]
    Q = lambda i: 1 + 1.0 * h[i[0]] ** 2 + 1e-4 * h[i[1]] ** 2
    r = run_protocol(_table(Q, [h, h]), [h, h], (1, 1), 1e-4, hs=[h, h], kinds=[DISC, DISC],
                     models=["power"] * 2, **V122)
    first = r["trace"][1]
    assert first[0] == 2 and first[1] == 1


def test_phase_two_can_coarsen_to_the_coarsest_value():
    # the second parameter has no effect: the recommendation loosens it all the way
    h = [0.4, 0.2, 0.1, 0.05, 0.025]
    t = [1e-2, 1e-4, 1e-6, 1e-8]
    Q = lambda i: 1 + h[i[0]] ** 2
    r = run_protocol(_table(Q, [h, t], cost=lambda i: float(4 ** i[0] * (1 + i[1]))),
                     [h, t], (2, 2), 1e-2, hs=[h, t], kinds=[DISC, NEST],
                     models=["power"] * 2, **V122)
    assert r["certified"] and r["proposal"][1] == 0


# ----------------------------------------------------------------------- job files
BASE = {"command": "echo {{n}}", "parameters": [{"name": "n", "values": [8, 16, 32]}],
        "observable": {"name": "u"}, "target": 1e-3}


@pytest.mark.parametrize("change, fragment", [
    ({"parameters": [{"name": "n", "knd": "nested-solver", "values": [8, 16]}]}, "'knd' is not a field"),
    ({"observable": {"name": "u", "regx": "u = (\\S+)"}}, "'regx' is not a field"),
    ({"cost": {"regx": "c=(\\S+)"}}, "'regx' is not a field"),
    ({"relative": "false"}, "'relative' must be true or false"),
    ({"keep_runs": "no"}, "'keep_runs' must be true or false"),
    ({"target": float("nan")}, "positive and finite"),
    ({"target": float("inf")}, "positive and finite"),
    ({"options": {"floor_probe": [{"n": 64}, {"m": 128}]}}, "give every parameter"),
    ({"failure": {"converged_if_regex": ["a", "b"]}}, "one pattern"),
    ({"render_files": "in.txt"}, "'render_files' must list"),
    ({"observable": {"name": "u", "column": 0}}, "counted from 1"),
    ({"observable": {"name": "u", "parse": "json"}}, "needs 'key'"),
    ({"parameters": [{"name": "run_id", "values": [1, 2]}], "command": "echo {{run_id}}"}, "reserved"),
    ({"env": ["A=1"]}, "'env' must map"),
    ({"options": "fast"}, "'options' must be a mapping"),
])
def test_job_file_refuses(change, fragment):
    with pytest.raises(JobError, match=None) as e:
        JobSpec.from_dict(dict(BASE, **change))
    assert fragment in str(e.value)


def test_column_zero_is_refused_when_reading():
    with pytest.raises(JobError, match="counted from 1"):
        extract_number({"column": 0}, "a 1 2 3")


def test_yaml_duplicate_keys_and_syntax_errors_are_job_errors(tmp_path):
    pytest.importorskip("yaml")
    p = tmp_path / "job.yaml"
    p.write_text("command: echo {{n}}\nparameters: [{name: n, values: [1, 2]}]\n"
                 "observable: u\ntarget: 1e-6\ntarget: 0.5\n")
    with pytest.raises(JobError, match="given twice"):
        load_job(str(p))
    p.write_text("command: [unclosed\n")
    with pytest.raises(JobError, match="not valid YAML"):
        load_job(str(p))


def test_listed_render_file_that_is_missing_fails_the_run(tmp_path):
    job = JobSpec.from_dict(dict(BASE, render_files=["inn.txt"]), base_dir=str(tmp_path))
    with pytest.raises(JobError, match="not in the run directory"):
        prepare_run_dir(job, str(tmp_path / "run"), {"n": 8})


def test_nan_observable_is_a_failed_run(tmp_path):
    cmd = f"{sys.executable} -c \"print('u = nan')\""
    job = JobSpec.from_dict(dict(BASE, command=cmd + " # {{n}}",
                                 observable={"name": "u", "regex": "u = (\\S+)"},
                                 workdir=str(tmp_path / "w")), base_dir=str(tmp_path))
    s = CommandSolver(job, log=None)
    with pytest.raises(RunFailed, match="not a finite number"):
        s(n=8)


def test_timed_out_run_is_retried_with_a_longer_limit(tmp_path):
    cmd = f"{sys.executable} -c \"import time; time.sleep(0.4); print('u = 1')\" # {{{{n}}}}"
    d = dict(BASE, command=cmd, observable={"name": "u", "regex": "u = (\\S+)"},
             workdir=str(tmp_path / "w"))
    short = JobSpec.from_dict(dict(d, timeout=0.1), base_dir=str(tmp_path))
    with pytest.raises(RunFailed, match="timed out"):
        CommandSolver(short, log=None)(n=8)
    longer = JobSpec.from_dict(dict(d, timeout=10), base_dir=str(tmp_path))
    out, _ = CommandSolver(longer, log=None)(n=8)
    assert out["u"] == 1.0
