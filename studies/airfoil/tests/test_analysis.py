"""The analyses on synthetic runs (analysis/synth.py): accounting, Stage 0, Stage 1, Stage 2.

The synthetic L/D has known errors per parameter, so the tests check that each
analysis finds what it must: the registered hardest design, the measured cost
model, the far-field and numerics rules, the reference and the admission of
targets, the grading (errbudget does not miss on this separable, power-law case),
setting (a), the optimum of the response surface and the payback.
"""
import contextlib
import io
import json
import math
import os
import shutil

import numpy as np
import pytest

import accounting
import replay
import stage0
import stage1
import stage2
import synth
from common import Store, engine, study, use_config

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(os.path.dirname(HERE), "config")


@pytest.fixture()
def lab(tmp_path):
    """A fresh store and configuration directory; restores the study's config dir afterwards."""
    cfg = tmp_path / "cfg"
    shutil.copytree(CONFIG, cfg)
    old = study.CONFIG_DIR
    use_config(str(cfg))
    yield tmp_path, str(tmp_path / "runs"), str(cfg)
    study.CONFIG_DIR = old


def run(main, *args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = main([str(x) for x in args])
    return code, buf.getvalue()


def set_result(runs, run_obj, **changes):
    p = os.path.join(runs, f"{run_obj.tag}.json")
    j = json.load(open(p))
    j.update(changes)
    json.dump(j, open(p, "w"))


# ---------------------------------------------------------------- accounting
def test_accounting_counts_and_resubmit_list(lab):
    tmp, runs, cfg = lab
    synth.write(runs, ["stage0"])
    s0 = study.sweep_stage0()
    os.remove(os.path.join(runs, f"{s0[3].tag}.json"))
    os.remove(os.path.join(runs, f"{s0[7].tag}.json"))
    with open(os.path.join(runs, f"{s0[7].tag}.error-1.json"), "w") as fh:
        json.dump(dict(stage="mesh", message="blockMesh exited with 1"), fh)
    for k in (1, 2, 3):
        with open(os.path.join(runs, f"{s0[9].tag}.error-{k}.json"), "w") as fh:
            json.dump(dict(stage="solver", message="x"), fh)
    os.remove(os.path.join(runs, f"{s0[9].tag}.json"))
    set_result(runs, s0[20], converged=False, stop_reason="stalled")
    a = accounting.account(Store(runs), "stage0")
    assert a["expected"] == 26 and a["missing"] == [3] and list(a["errors"]) == [7] and a["gave_up"] == [9]
    assert a["stalled"] == [20] and a["results"] == 23 and not a["complete"]
    assert a["resubmit"] == "3,7"
    sacct = tmp / "sacct.txt"
    sacct.write_text("JobID|State\n100_0|COMPLETED\n100_0.batch|COMPLETED\n100_1|TIMEOUT\n100_2|OUT_OF_MEMORY\n")
    assert accounting.sacct_states(str(sacct)) == {"COMPLETED": 1, "TIMEOUT": 1, "OUT_OF_MEMORY": 1}


# ---------------------------------------------------------------- Stage 0
def test_stage0_cost_model_hardest_and_config(lab):
    tmp, runs, cfg = lab
    synth.write(runs, ["stage0"])
    code, out = run(stage0.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res"), "--write-config")
    assert code == 0, out
    s = json.load(open(tmp / "res" / "stage0_summary.json"))
    assert s["cost_model"]["t_cell_iter"]["median"] == pytest.approx(synth.T_TRUE, rel=1e-6)
    assert s["cost_model"]["iter_growth_per_level"] == pytest.approx(study.ITER_GROWTH, rel=0.02)
    assert s["cost_model"]["parallel"]["speedup"] == pytest.approx(3.4, rel=1e-6)
    # (6, 8) and (6, 16) tie within 10 %; the tie goes to the thickest
    assert s["hardest"] == list(study.CORNERS[-1]) and s["radius_factor"] == 1 and s["numerics_check_triggered"] is False
    assert json.load(open(os.path.join(cfg, "stage1.json")))["hardest"] == list(study.CORNERS[-1])
    assert study.hardest_design() == study.CORNERS[-1]
    assert 50 < s["predicted"]["median"]["stage1"] < 200


def test_stage0_rules_failure_far_field_fallback(lab):
    tmp, runs, cfg = lab
    synth.write(runs, ["stage0"], far=0.25)               # far-field error ten times larger
    g = stage0.runs_by_group()
    r181 = next(r for r in g["C"] if (r.m, r.t) == (0.0, 0.08) and r.n_surf == 181)
    set_result(runs, r181, converged=False, stop_reason="max_iter")
    store = Store(runs)
    design, rows = stage0.hardest(store, g["C"])
    assert design == [0.0, 0.08]                          # an unconverged corner ranks first
    factor, frows = stage0.far_field(store, g["D"])
    # A5: the radius list never changes; the far-field error at R = 5000 is reported, and with
    # the synthetic 1/R error it is exactly the model's term at R = 5000
    assert factor == 1 and all(r["above_eps10"] for r in frows)
    for row in frows:
        q = row["ld"]
        assert row["delta_R"] == pytest.approx((q[5000.0] - q[500.0]) / 9, rel=1e-12)
    fb, info = stage0.numerics_check(store, study.sweep_stage0())
    assert not fb
    for r in g["D"][:3]:
        set_result(runs, r, converged=False, stop_reason="stalled")
    fb, info = stage0.numerics_check(Store(runs), study.sweep_stage0())
    assert fb and len(info["unconverged_at_1e6"]) == 3
    # A3: the check triggering stops the study; no configuration is written, no fallback
    code, out = run(stage0.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res3"), "--write-config")
    assert code == 1 and "STOP (PROTOCOL.md 3.3, A3)" in out and "NOT writing config" in out
    nj = os.path.join(cfg, "numerics.json")
    assert not os.path.exists(nj) or json.load(open(nj)).get("numerics") != "fallback"
    # an incomplete Stage 0 never writes the configuration
    os.remove(os.path.join(runs, f"{g['A'][0].tag}.json"))
    code, out = run(stage0.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res"), "--write-config")
    assert code == 1 and "NOT writing config" in out
    assert json.load(open(os.path.join(cfg, "stage1.json")))["hardest"] is None


def test_radius_list_never_changes_after_a5(lab):
    """A5 withdrew the radius amendment of PROTOCOL.md 3.2: a factor in stage1.json is ignored."""
    tmp, runs, cfg = lab
    json.dump(dict(hardest=list(study.CORNERS[2]), radius_factor=10), open(os.path.join(cfg, "stage1.json"), "w"))
    assert study.ladders_now()[2] == [50.0, 500.0, 5000.0]
    assert study.reference_setting()["radius"] == 5000.0
    s1 = study.sweep_stage1()
    assert len(s1) == 540 and {r.radius for r in s1} == {50.0, 500.0, 5000.0}
    ref = study.sweep_stage1ref()
    assert all(r.settle for r in ref) and [r.radius for r in ref] == [5000.0, 5000.0, 5000.0, 500.0, 5000.0, 500.0]
    assert ref[-1].n_surf == 91 and ref[-1].res_ctrl == 1e-6          # the far-field run


# ---------------------------------------------------------------- the replay copy
def test_replay_copy_matches_the_article_baselines():
    """seq/coord/ref_check are verbatim copies of the article's paper/baselines.py."""
    import inspect
    src = os.path.join(os.path.expanduser("~/OpenResearch/knobs-paper/knobs-paper/paper"), "baselines.py")
    if not os.path.exists(src):
        pytest.skip("article repository not on this machine")
    text = open(src).read()
    for fn in (replay.seq, replay.coord, replay.ref_check):
        assert inspect.getsource(fn) in text, fn.__name__


# ---------------------------------------------------------------- Stage 1
def _stage1_data(runs, cfg, scale):
    json.dump(dict(hardest=list(study.CORNERS[2]), radius_factor=1), open(os.path.join(cfg, "stage1.json"), "w"))
    synth.write(runs, ["stage1", "stage1ref"], scale=scale)


def test_stage1_admitted_grading_and_settings(lab):
    tmp, runs, cfg = lab
    _stage1_data(runs, cfg, scale=0.1)                    # errors small enough for 1 % to be admitted
    code, out = run(stage1.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res"),
                    "--starts-every", "25", "--write-settings")
    assert code == 0, out
    s = json.load(open(tmp / "res" / "stage1_summary.json"))
    ref = s["reference"]
    # order 2 by construction; the +-2e-6 iterative error of the synthetic runs blurs it a little
    assert ref["u"] < 1e-3 and ref["mesh_order"] == pytest.approx(2.0, abs=0.3)
    # A5: the synthetic far-field error is exactly 1/R, so L/D(d) + Delta_R removes it: the
    # reference value is the model's value without its far-field term (to the iterative noise)
    d = study.stage1ref_runs(study.hardest_design())[0]
    no_far = synth.model(d, scale=0.1, far=0.0)
    assert ref["reference_value"] == pytest.approx(no_far, rel=1e-4)
    # it corrected by the model's far-field term at R = 5000 (2.5e-5 of L/D at this scale), to
    # within the synthetic iterative error of the two n_surf 91 runs (+-2e-5 each, over 9)
    assert ref["reference_value"] - ref["values"]["d"] == pytest.approx(
        no_far - synth.model(d, scale=0.1), abs=2 * 2e-5 * abs(no_far) / 9 * 1.01)
    ff = ref["far_field"]
    assert ff["delta"] == pytest.approx((ff["ld_5000"] - ff["ld_500"]) / 9, rel=1e-12)
    assert ref["u_far_field"] == pytest.approx(0.2 * abs(ff["delta"]) / abs(ref["reference_value"]), rel=1e-12)
    assert ref["u_adm"] >= ref["u"] + ref["u_far_field"] - 1e-15
    assert [t for t, adm in s["targets"] if adm][0] == 0.01
    g = s["grades"]["0.01"]["rows"]
    assert g["program"]["missed"] == 0 and g["program"]["reached"] > 0
    assert g["ref"]["reached"] + g["ref"]["missed"] == s["grades"]["0.01"]["n_starts"]
    uc = s["user_case"]["0.01"]
    assert uc["certified"] and uc["error_vs_reference"] <= 0.01 and uc["K"] > 1 and uc["runs_used"] > 5
    st = json.load(open(os.path.join(cfg, "stage2_settings.json")))["settings"]
    assert sorted(st) == ["a", "b", "c", "d"]
    assert st["b"] == dict(study.default_start(), res_ctrl=1e-5, gamg_tol=1e-6)
    assert st["c"] == study.cautious_setting() and st["d"] == study.reference_setting()
    assert len(study.sweep_stage2()) == 125                # A5: 25 far-field runs with (d)


def test_stage1_not_admitted_falls_back_to_three_percent(lab):
    tmp, runs, cfg = lab
    _stage1_data(runs, cfg, scale=1.0)
    code, out = run(stage1.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res"), "--starts-every", "40")
    assert code == 0, out
    s = json.load(open(tmp / "res" / "stage1_summary.json"))
    adm = dict((t, a) for t, a in s["targets"])
    assert adm[0.01] is False and adm[0.03] is True
    assert s["grades"]["0.03"]["rows"]["program"]["missed"] == 0
    assert "not admitted" in out


def test_unconverged_runs_are_failures_for_every_procedure(lab):
    tmp, runs, cfg = lab
    _stage1_data(runs, cfg, scale=0.1)
    T = stage1.tensor_of(Store(runs), study.sweep_stage1())
    bad = [i for i in T if not T[i][2]["converged"]]
    assert len(bad) == 18                                  # 9 stalled + 9 diverged by construction
    eb = engine()
    lad = study.ladders_now()
    res = replay.grade(T, lad, study.hs_of(lad), study.KINDS, 80.0, 0.8, eb.run_protocol, dict(eb.V121),
                       starts=[study.GRID[0], study.GRID[100]])
    for p, o in res["procs"].items():
        for rec in o["records"]:
            if p == "program":
                assert rec["returned"] is None or rec["returned"] not in bad
            elif rec["returned"] in bad:
                assert rec["error"] is None and any(m[1] == rec["returned"] for m in o["misses"])


# ---------------------------------------------------------------- Stage 2
def test_response_surface_fit_and_box_optimum():
    designs = study.DESIGNS
    c_true = np.array([90.0, 2.0, -1.0, -3.0, 0.5, -2.0])
    vals = [float(stage2.surface(c_true, m, t)) for m, t in designs]
    c, rms, r2 = stage2.fit(designs, vals)
    assert np.allclose(c, c_true) and rms < 1e-10 and r2 == pytest.approx(1.0)
    m, t, f = stage2.optimum(c_true)
    x, y = np.linalg.solve([[-6.0, 0.5], [0.5, -4.0]], [-2.0, 1.0])   # stationary point of the scaled quadratic
    assert m == pytest.approx(stage2.M0 + stage2.MS * x, abs=2e-4) and t == pytest.approx(stage2.T0 + stage2.TS * y, abs=2e-4)
    m, t, f = stage2.optimum(np.array([0.0, 1.0, 1.0, 0.0, 0.0, 0.0]))          # increasing: corner
    assert (m, t) == (max(study.CAMBERS), max(study.THICKNESSES))


def test_stage2_optimum_costs_and_payback(lab):
    tmp, runs, cfg = lab
    _stage1_data(runs, cfg, scale=0.1)
    run(stage1.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res"), "--starts-every", "25",
        "--write-settings")
    n = synth.write(runs, ["stage2"], scale=0.1)
    assert n < 125                                        # the hardest design's runs are reused
    code, out = run(stage2.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res"),
                    "--stage1", str(tmp / "res" / "stage1_summary.json"))
    assert code == 0, out
    s = json.load(open(tmp / "res" / "stage2_summary.json"))
    d = s["per_setting"]["d"]["opt"]
    assert d["m"] == pytest.approx(0.04, abs=1e-3) and d["t"] == pytest.approx(0.10, abs=1e-3)
    a = s["per_setting"]["a"]["against_reference"]
    assert a["loss"] < 1e-3 and a["max_error"] < 0.01
    pb = s["payback"]["c"]
    assert pb["n_star"] == pytest.approx(pb["K_prime"] / (pb["mu"] - 1.0))
    assert s["search"]["runs"] > 5


# ---------------------------------------------------------------- the live job file
def test_live_job_file_returns_the_replay_setting(lab):
    """errbudget_job.yaml through live_tune.py (the command interface, driver.py --print
    on the stored runs) must return the setting the Stage-1 replay returns."""
    import subprocess
    import sys
    tmp, runs, cfg = lab
    _stage1_data(runs, cfg, scale=0.1)
    run(stage1.main, "--runs", runs, "--config", cfg, "--out", str(tmp / "res"), "--starts-every", "60")
    want = json.load(open(tmp / "res" / "stage1_summary.json"))["user_case"]["0.01"]
    root = os.path.dirname(HERE)
    env = dict(os.environ, AIRFOIL_WORK=str(tmp), AIRFOIL_CONFIG=cfg, PYTHONDONTWRITEBYTECODE="1",
               PATH=f"{os.path.dirname(sys.executable)}:{os.environ['PATH']}")
    r = subprocess.run([sys.executable, os.path.join(root, "live_tune.py"), os.path.join(root, "errbudget_job.yaml"),
                        "--out", str(tmp / "live")], env=env, capture_output=True, text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    got = json.load(open(tmp / "live" / "errbudget_result.json"))
    assert got["certified"] == want["certified"]
    assert {k: float(v) for k, v in got["settings"].items()} == {k: float(v) for k, v in want["setting"].items()}
