"""The driver end to end, against the fake OpenFOAM of tests/fakefoam (no OpenFOAM needed).

Covers: a converged run and its JSON, idempotence, a stalled run, a diverged run,
a set-up error retried and given up, an interrupt (SIGUSR1) and resume, a run on two
cores, the single-run interface of the errbudget job file, --list/--count/--groups,
and the dry run of every stage.
"""
import json
import os
import signal
import subprocess
import sys
import time

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DRIVER = os.path.join(ROOT, "driver.py")
FAKE = os.path.join(HERE, "fakefoam", "fakefoam.py")
TOOLS = ("checkMesh", "simpleFoam", "mapFields", "decomposePar", "reconstructPar", "mpirun")
STAGE2 = {"a": {"n_surf": 91, "dy1": 2e-5, "radius": 500.0, "res_ctrl": 1e-4, "gamg_tol": 1e-7},
          "b": {"n_surf": 64, "dy1": 2e-5, "radius": 500.0, "res_ctrl": 1e-5, "gamg_tol": 1e-6},
          "c": {"n_surf": 256, "dy1": 1e-5, "radius": 5000.0, "res_ctrl": 1e-6, "gamg_tol": 1e-8},
          "d": {"n_surf": 362, "dy1": 5e-6, "radius": 50000.0, "res_ctrl": 1e-7, "gamg_tol": 1e-9}}


@pytest.fixture()
def env(tmp_path):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    os.chmod(FAKE, 0o755)
    for t in TOOLS:
        (bindir / t).symlink_to(FAKE)
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    (cfg / "stage1.json").write_text(json.dumps({"hardest": [0.06, 0.08], "radius_factor": 1}))
    (cfg / "stage2_settings.json").write_text(json.dumps({"settings": STAGE2}))
    (cfg / "numerics.json").write_text(json.dumps({"numerics": "primary"}))
    e = dict(os.environ)
    e.update(PATH=f"{bindir}:{os.path.dirname(sys.executable)}:{e['PATH']}",
             AIRFOIL_WORK=str(tmp_path / "work"), AIRFOIL_CONFIG=str(cfg), AIRFOIL_POLL_S="0.05",
             AIRFOIL_STALL_ITER="60", PYTHONDONTWRITEBYTECODE="1", SLURM_CPUS_PER_TASK="2")
    for k in ("FAKEFOAM_STALL", "FAKEFOAM_CRASH_AT", "FAKEFOAM_SETUP_FAIL", "FAKEFOAM_SLEEP"):
        e.pop(k, None)
    return e


def drive(env, *args, **extra):
    e = dict(env, **extra)
    return subprocess.run([sys.executable, DRIVER, *args], env=e, capture_output=True, text=True)


def result(env, tag):
    return json.load(open(os.path.join(env["AIRFOIL_WORK"], "runs", f"{tag}.json")))


SMALL = ["m=0.03", "t=0.12", "n_surf=64", "dy1=2e-5", "radius=50", "gamg_tol=1e-7"]


def tag_of(env, *kv):
    r = drive(env, "--run-setting", *kv, "--list")
    return r.stdout.split()[1]


def test_converged_run_and_idempotence(env):
    r = drive(env, "--sweep", "check", "--index", "0")
    assert r.returncode == 0, r.stdout + r.stderr
    tag = "m2.0_t10_n64_y2e-05_R50_rc1e-03_gt1e-07"
    j = result(env, tag)
    assert j["status"] == "ok" and j["converged"] and j["stop_reason"] == "converged"
    assert 0 < j["ld"] < 200 and j["cl"] > 0 and j["cd"] > 0
    assert j["cost"]["cells_iters"] == j["cells"] * j["iterations"]
    assert j["cells"] == j["mesh"]["generator"]["cells"]
    assert j["yplus"]["n_faces"] == 128 and j["yplus"]["max"] > j["yplus"]["min"] > 0
    assert j["first_cell"]["median"] > 0 and j["openfoam"]["version"] == "v2606"
    assert j["sweeps"] == [["check", 0]]
    assert os.path.exists(os.path.join(env["AIRFOIL_WORK"], "runs", f"{tag}.log.txt"))
    assert not os.path.exists(os.path.join(env["AIRFOIL_WORK"], "cases", tag))
    r2 = drive(env, "--sweep", "check", "--index", "0")
    assert r2.returncode == 0 and "done already" in r2.stdout
    st = drive(env, "--sweep", "check", "--status")
    assert "done 1" in st.stdout and "missing 1" in st.stdout


def test_stalled_run(env):
    r = drive(env, "--run-setting", *SMALL, "res_ctrl=1e-6", FAKEFOAM_STALL="1")
    assert r.returncode == 0, r.stdout + r.stderr
    j = result(env, tag_of(env, *SMALL, "res_ctrl=1e-6"))
    assert j["status"] == "ok" and not j["converged"] and j["stop_reason"] == "stalled"
    assert j["iterations"] < 2000


def test_diverged_run_is_a_failed_result(env):
    r = drive(env, "--run-setting", *SMALL, "res_ctrl=1e-5", FAKEFOAM_CRASH_AT="20")
    assert r.returncode == 0, r.stdout + r.stderr
    j = result(env, tag_of(env, *SMALL, "res_ctrl=1e-5"))
    assert j["status"] == "failed" and not j["converged"] and j["stop_reason"] == "diverged"
    p = drive(env, "--run-setting", *SMALL, "res_ctrl=1e-5", "--print")
    assert p.returncode == 1 and "RUN FAILED" in p.stdout


def test_setup_error_retried_then_given_up(env):
    for k in range(3):
        r = drive(env, "--run-setting", *SMALL, "res_ctrl=1e-4", FAKEFOAM_SETUP_FAIL="1")
        assert r.returncode == 2, r.stdout + r.stderr
    errs = [f for f in os.listdir(os.path.join(env["AIRFOIL_WORK"], "runs")) if ".error-" in f]
    assert len(errs) == 3
    r = drive(env, "--run-setting", *SMALL, "res_ctrl=1e-4")
    assert r.returncode == 4 and "GAVE UP" in r.stdout


def test_interrupt_and_resume(env):
    kv = [*SMALL, "res_ctrl=1e-6"]
    e = dict(env, FAKEFOAM_SLEEP="0.01")
    p = subprocess.Popen([sys.executable, DRIVER, "--run-setting", *kv], env=e,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    time.sleep(3.0)
    p.send_signal(signal.SIGUSR1)
    out, _ = p.communicate(timeout=120)
    assert p.returncode == 3, out
    tag = tag_of(env, *kv)
    state = json.load(open(os.path.join(env["AIRFOIL_WORK"], "cases", tag, "STATE.json")))
    assert state["segments"] == 1 and state["resume_from"] > 0
    r = drive(env, "--run-setting", *kv)
    assert r.returncode == 0 and "resuming" in r.stdout, r.stdout + r.stderr
    j = result(env, tag)
    assert j["converged"] and j["segments"] == 2
    # the fake converges where 0.2 * 0.985**it < 1e-6: the resumed run must end there too
    import math
    assert j["iterations"] == math.ceil(math.log(1e-6 / 0.2) / math.log(0.985) + 1e-9)


def test_parallel_run(env):
    r = drive(env, "--run-setting", *SMALL, "res_ctrl=1e-3", "np=2")
    assert r.returncode == 0, r.stdout + r.stderr
    j = result(env, tag_of(env, *SMALL, "res_ctrl=1e-3", "np=2"))
    assert j["tag"].endswith("_np2")
    assert j["np"] == 2 and j["cost"]["core_s"] == pytest.approx(2 * j["cost"]["cpu_s"])
    assert j["yplus"]["n_faces"] == 128


def test_print_interface_for_errbudget(env):
    kv = [*SMALL, "res_ctrl=1e-3"]
    r = drive(env, "--run-setting", *kv, "--print")
    assert r.returncode == 0, r.stdout + r.stderr
    vals = dict(line.split(" = ") for line in r.stdout.splitlines() if " = " in line)
    assert float(vals["L/D"]) > 0 and float(vals["cost"]) > 0 and "converged = yes" in r.stdout
    r = drive(env, "--run-setting", "n_surf=64", "dy1=2e-5", "radius=50", "res_ctrl=1e-3",
              "gamg_tol=1e-7", "--design-from-config", "--print")
    assert r.returncode == 0 and "L/D = " in r.stdout


def test_list_count_groups(env):
    for sweep, n in (("check", 2), ("meshcheck", 6), ("stage0", 26), ("stage1", 540), ("stage1ref", 6),
                     ("stage2", 125)):
        r = drive(env, "--sweep", sweep, "--count")
        assert r.returncode == 0 and int(r.stdout) == n
        lst = drive(env, "--sweep", sweep, "--list", "--json")
        rows = [json.loads(x) for x in lst.stdout.splitlines()]
        assert len(rows) == n and len({x["tag"] for x in rows}) == n
    g = drive(env, "--sweep", "stage1", "--groups").stdout.split("\n")
    assert g[0].split() == ["03:00:00", "2G", "1", "0-323"]
    assert g[1].split() == ["06:00:00", "3G", "1", "324-431"]
    assert g[2].split() == ["12:00:00", "4G", "1", "432-539"]


@pytest.mark.parametrize("sweep,index", [("check", 1), ("stage0", 25), ("stage1", 539),
                                         ("stage1ref", 4), ("stage2", 99)])
def test_dry_run_writes_complete_case(env, tmp_path, sweep, index):
    case = tmp_path / f"dry_{sweep}"
    r = drive(env, "--sweep", sweep, "--index", str(index), "--dry-run", "--case-dir", str(case))
    assert r.returncode == 0, r.stdout + r.stderr
    for f in ("system/controlDict", "system/fvSchemes", "system/fvSolution",
              "system/decomposeParDict", "constant/transportProperties", "constant/turbulenceProperties",
              "constant/polyMesh/points", "constant/polyMesh/faces", "constant/polyMesh/owner",
              "constant/polyMesh/neighbour", "constant/polyMesh/boundary",
              "0/U", "0/p", "0/k", "0/omega", "0/nut"):
        assert (case / f).stat().st_size > 100, f
    assert not (case / "system/blockMeshDict").exists()


def test_result_of_an_older_mesh_is_set_aside_and_rerun(env):
    r = drive(env, "--sweep", "check", "--index", "0")
    assert r.returncode == 0
    tag = "m2.0_t10_n64_y2e-05_R50_rc1e-03_gt1e-07"
    p = os.path.join(env["AIRFOIL_WORK"], "runs", f"{tag}.json")
    j = json.load(open(p))
    j.pop("mesh_version")                                   # as written by the first generator
    json.dump(j, open(p, "w"))
    st = drive(env, "--sweep", "check", "--status")
    assert "stale 1" in st.stdout
    assert "0-1" in drive(env, "--sweep", "check", "--groups").stdout   # resubmitted by submit.sh
    r = drive(env, "--sweep", "check", "--index", "0")
    assert r.returncode == 0 and "older mesh" in r.stdout
    assert os.path.exists(os.path.join(env["AIRFOIL_WORK"], "runs", "stale", f"{tag}.json"))
    assert result(env, tag)["mesh_version"] >= 2


def test_result_of_an_older_set_up_is_set_aside_and_rerun(env):
    """A3: results made with the registered free stream and numerics (no setup_version)
    are set aside like results of an older mesh."""
    assert drive(env, "--sweep", "check", "--index", "0").returncode == 0
    tag = "m2.0_t10_n64_y2e-05_R50_rc1e-03_gt1e-07"
    p = os.path.join(env["AIRFOIL_WORK"], "runs", f"{tag}.json")
    j = json.load(open(p))
    assert j["setup_version"] == 2
    j.pop("setup_version")
    json.dump(j, open(p, "w"))
    assert "stale 1" in drive(env, "--sweep", "check", "--status").stdout
    r = drive(env, "--sweep", "check", "--index", "0")
    assert r.returncode == 0 and "set aside" in r.stdout
    assert result(env, tag)["setup_version"] == 2


def test_stop_followed_by_an_io_error_is_the_stop_not_a_divergence(env):
    """On Hamilton the stall watchdog's stop met a fatal I/O error (controlDict re-read while
    being written) and the run was filed as diverged.  Now the stop is what is recorded."""
    kv = [*SMALL, "res_ctrl=1e-6"]
    r = drive(env, "--run-setting", *kv, FAKEFOAM_STALL="1", FAKEFOAM_IO_FAIL_ON_STOP="1")
    assert r.returncode == 0, r.stdout + r.stderr
    j = result(env, tag_of(env, *kv))
    assert j["status"] == "ok" and j["stop_reason"] == "stalled" and not j["converged"]


def test_per_run_numerics_and_diag_sweep(env):
    rows = [json.loads(x) for x in drive(env, "--sweep", "diag", "--list", "--json").stdout.splitlines()]
    assert len(rows) == 5 and rows[0]["tag"].endswith("_fb") and not rows[1]["tag"].endswith("_fb")
    case = os.path.join(env["AIRFOIL_WORK"], "dry_fb")
    r = drive(env, "--sweep", "diag", "--index", "0", "--dry-run", "--case-dir", case)
    assert r.returncode == 0
    fv = open(os.path.join(case, "system", "fvSolution")).read()
    fs = open(os.path.join(case, "system", "fvSchemes")).read()
    # the registered fallback (withdrawn by A3, kept for the diag sweep)
    assert "consistent      no;" in fv and "p               0.3;" in fv and "bounded Gauss upwind" in fs


def test_failed_check_mesh_is_an_error_not_a_result(env):
    r = drive(env, "--run-setting", *SMALL, "res_ctrl=1e-3", FAKEFOAM_CHECKMESH_FAIL="1")
    assert r.returncode == 2, r.stdout + r.stderr
    errs = [f for f in os.listdir(os.path.join(env["AIRFOIL_WORK"], "runs")) if ".error-" in f]
    e = json.load(open(os.path.join(env["AIRFOIL_WORK"], "runs", errs[0])))
    assert e["stage"] == "mesh" and "face pyramids" in e["message"]
    assert not [f for f in os.listdir(os.path.join(env["AIRFOIL_WORK"], "runs")) if f.endswith(".json") and ".error-" not in f]


def test_mesh_only_check(env):
    r = drive(env, "--sweep", "meshcheck", "--index", "0", "--mesh-only")
    assert r.returncode == 0 and "MESH PASSED" in r.stdout, r.stdout + r.stderr
    d = os.path.join(env["AIRFOIL_WORK"], "meshcheck")
    js = [f for f in os.listdir(d) if f.endswith(".json")]
    m = json.load(open(os.path.join(d, js[0])))
    assert m["passed"] and m["checkMesh"]["cells"] == m["cells"] and m["checkMesh"]["failures"] == []
    assert not os.listdir(os.path.join(env["AIRFOIL_WORK"], "cases"))      # cleaned up
    r = drive(env, "--sweep", "meshcheck", "--index", "1", "--mesh-only", FAKEFOAM_CHECKMESH_FAIL="1")
    assert r.returncode == 2 and "MESH FAILED" in r.stdout


def lock_file(env, tag):
    return os.path.join(env["AIRFOIL_WORK"], "runs", f"{tag}.running")


def test_second_copy_of_a_running_run_does_nothing(env):
    """10 Nov 2026: the diag sweep, submitted again while two of its runs still ran, started
    a second copy of each in the same case directory; the first copy to finish deleted the
    case under the other.  Now the second copy finds the lock and leaves."""
    tag = "m2.0_t10_n64_y2e-05_R50_rc1e-03_gt1e-07"
    e = dict(env, FAKEFOAM_SLEEP="0.01")
    a = subprocess.Popen([sys.executable, DRIVER, "--sweep", "check", "--index", "0"], env=e,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for _ in range(400):
        if os.path.exists(os.path.join(env["AIRFOIL_WORK"], "cases", tag, "log.simpleFoam.0")):
            break
        time.sleep(0.05)
    b = drive(env, "--sweep", "check", "--index", "0")
    assert b.returncode == 5 and "RUNNING ELSEWHERE" in b.stdout, b.stdout + b.stderr
    assert "running 1" in drive(env, "--sweep", "check", "--status").stdout
    assert drive(env, "--sweep", "check", "--groups").stdout.split()[-1] == "1"   # only index 1 to submit
    out, _ = a.communicate(timeout=300)
    assert a.returncode == 0, out
    assert result(env, tag)["status"] == "ok" and not os.path.exists(lock_file(env, tag))


def test_lock_of_a_dead_process_is_taken_over(env):
    tag = "m2.0_t10_n64_y2e-05_R50_rc1e-03_gt1e-07"
    os.makedirs(os.path.join(env["AIRFOIL_WORK"], "runs"))
    dead = subprocess.Popen([sys.executable, "-c", "pass"])
    dead.wait()
    import socket
    json.dump(dict(host=socket.gethostname(), pid=dead.pid, slurm={}, time="earlier"),
              open(lock_file(env, tag), "w"))
    assert "running" not in drive(env, "--sweep", "check", "--status").stdout
    r = drive(env, "--sweep", "check", "--index", "0")
    assert r.returncode == 0, r.stdout + r.stderr
    assert not os.path.exists(lock_file(env, tag))


def test_numerics_sweep_keeps_fields_and_wall_diagnostics(env):
    import numpy as np
    rows = [json.loads(x) for x in drive(env, "--sweep", "numerics", "--list", "--json").stdout.splitlines()]
    assert len(rows) == 33 and len({x["tag"] for x in rows}) == 33
    i = next(x["index"] for x in rows if x["tag"].endswith("_nvrobust") and "_n64_" in x["tag"])
    tag = rows[i]["tag"]
    r = drive(env, "--sweep", "numerics", "--index", str(i))
    assert r.returncode == 0, r.stdout + r.stderr
    j = result(env, tag)
    assert j["errors"] == [], j["errors"]
    f = j["fields"]
    assert set(f["wall"]) == {"upper", "lower"} and len(f["wall"]["upper"]["x"]) == 64
    assert len(f["wall"]["lower"]["cp"]) == 64 and f["k_max"]["value"] > 0
    assert f["change"]["max_per_iter"] > 0 and abs(sum(f["change"]["share"].values()) - 1) < 1e-9
    assert f["change"]["between"][1] - f["change"]["between"][0] > 0
    assert {"it", "Ux", "Uy", "p", "k", "omega"} <= set(j["residual_fields"])
    z = np.load(os.path.join(env["AIRFOIL_WORK"], "runs", "fields", f"{tag}.npz"))
    assert z["U"].shape == (j["cells"], 2) and z["G"].shape[2] == 2 and z["dU"].shape == (j["cells"],)


def test_numerics_variants_write_their_dictionaries(env, tmp_path):
    rows = [json.loads(x) for x in drive(env, "--sweep", "numerics", "--list", "--json").stdout.splitlines()]
    sa = next(x["index"] for x in rows if x["tag"].endswith("_nvsa"))
    nonorth = next(x["index"] for x in rows if x["tag"].endswith("_nvrobust_nonorth"))
    case = tmp_path / "sa"
    assert drive(env, "--sweep", "numerics", "--index", str(sa), "--dry-run", "--case-dir", str(case)).returncode == 0
    assert "SpalartAllmaras" in (case / "constant/turbulenceProperties").read_text()
    assert "div(phi,nuTilda) bounded Gauss upwind" in (case / "system/fvSchemes").read_text()
    assert '"(U|nuTilda)"' in (case / "system/fvSolution").read_text()
    assert (case / "0/nuTilda").exists()
    ctl = (case / "system/controlDict").read_text()
    assert "writeInterval   100;" in ctl and "purgeWrite      2;" in ctl
    case = tmp_path / "nonorth"
    assert drive(env, "--sweep", "numerics", "--index", str(nonorth), "--dry-run", "--case-dir", str(case)).returncode == 0
    assert "nNonOrthogonalCorrectors 1;" in (case / "system/fvSolution").read_text()
    assert not (case / "0/nuTilda").exists()
    assert "grad(U)         cellLimited Gauss linear 1;" in (case / "system/fvSchemes").read_text()


def test_sequenced_start_of_a_fine_run(env):
    """A4: a run finer than the Stage-1 lists first solves the same design on the coarse
    sequence mesh of its radius, maps it, and carries its cost (the threshold is lowered for
    the test so that n_surf 128 counts as fine)."""
    e = dict(env, AIRFOIL_SEQUENCE_ABOVE="91")
    kv = ["m=0.04", "t=0.12", "n_surf=128", "dy1=2e-5", "radius=50", "res_ctrl=1e-6", "gamg_tol=1e-8"]
    r = drive(e, "--run-setting", *kv)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "first run m4.0_t12_n91_y2e-05_R50_rc1e-06_gt1e-08_st" in r.stdout
    j = result(env, "m4.0_t12_n128_y2e-05_R50_rc1e-06_gt1e-08")
    s = j["sequence"]
    assert s["tag"] == "m4.0_t12_n91_y2e-05_R50_rc1e-06_gt1e-08_st" and s["converged"] and s["iterations"] > 0
    assert s["ld"] > 0                                    # A5: the first run's L/D, for the correction
    assert j["cost"]["cells_iters"] == j["cells"] * j["iterations"] + s["cells"] * s["iterations"]
    assert j["cost"]["cpu_s"] > s["cpu_s"] > 0
    log = open(os.path.join(env["AIRFOIL_WORK"], "runs", f"{j['tag']}.log.txt")).read()
    assert "log.sequence.simpleFoam" in log and "log.mapFields" in log
    assert not os.listdir(os.path.join(env["AIRFOIL_WORK"], "cases"))
    # a coarse run is not sequenced
    r = drive(e, "--run-setting", "m=0.04", "t=0.12", "n_surf=64", "dy1=2e-5", "radius=50", "res_ctrl=1e-3", "gamg_tol=1e-7")
    assert r.returncode == 0 and "first run" not in r.stdout
    assert "sequence" not in result(env, "m4.0_t12_n64_y2e-05_R50_rc1e-03_gt1e-07")


def test_diverging_first_run_is_an_error(env):
    e = dict(env, AIRFOIL_SEQUENCE_ABOVE="91", FAKEFOAM_CRASH_AT="30")
    kv = ["m=0.04", "t=0.12", "n_surf=128", "dy1=2e-5", "radius=50", "res_ctrl=1e-6", "gamg_tol=1e-8"]
    r = drive(e, "--run-setting", *kv)
    assert r.returncode == 2, r.stdout + r.stderr
    errs = [f for f in os.listdir(os.path.join(env["AIRFOIL_WORK"], "runs")) if ".error-" in f]
    msg = json.load(open(os.path.join(env["AIRFOIL_WORK"], "runs", errs[0])))["message"]
    assert "first run" in msg and "diverged" in msg


def test_settle_runs_ignore_the_stall_watchdog_and_keep_their_history(env):
    """The settling measurement: no residual stop (residualControl 1e-12), no stall watchdog,
    the whole force history in runs/history/<tag>.csv.  A small stand-in mesh replaces the
    reference mesh, and the fake residuals plateau (a normal run would stall after 60)."""
    e = dict(env, AIRFOIL_SETTLE_TEST_MESH="64,4e-5", AIRFOIL_SEQUENCE_ABOVE="1000", FAKEFOAM_STALL="1")
    rows = [json.loads(x) for x in drive(e, "--sweep", "settle", "--list", "--json").stdout.splitlines()]
    assert len(rows) == 2 and all("rc1e-12" in x["tag"] for x in rows) and rows[0]["time"] == "2-12:00:00"
    r = drive(e, "--sweep", "settle", "--index", "1")
    assert r.returncode == 0, r.stdout + r.stderr
    j = result(env, rows[1]["tag"])
    assert j["stop_reason"] == "max_iter" and j["iterations"] == 40000 and not j["converged"]
    h = os.path.join(env["AIRFOIL_WORK"], "runs", j["history_file"])
    lines = open(h).read().splitlines()
    assert lines[0] == "iteration,cl,cd,cm" and len(lines) == 40001


def test_settled_run_is_stopped_by_the_driver(env, tmp_path):
    """A5: a settled run has OpenFOAM's residual stop off (residualControl 1e-12) and is
    stopped by the driver once its largest residual is at its target and L/D has settled
    over the window (lowered to 200 iterations for the fake solver)."""
    rows = [json.loads(x) for x in drive(env, "--sweep", "stage1ref", "--list", "--json").stdout.splitlines()]
    ff = rows[-1]
    assert ff["tag"].endswith("_R500_rc1e-06_gt1e-08_st") and "_n91_" in ff["tag"]
    case = tmp_path / "ff"
    assert drive(env, "--sweep", "stage1ref", "--index", str(ff["index"]), "--dry-run", "--case-dir", str(case)).returncode == 0
    assert "p               1e-12;" in (case / "system/fvSolution").read_text()
    r = drive(dict(env, AIRFOIL_SETTLE_WINDOW="200"), "--sweep", "stage1ref", "--index", str(ff["index"]))
    assert r.returncode == 0, r.stdout + r.stderr
    j = result(env, ff["tag"])
    assert j["stop_reason"] == "settled" and j["converged"] and j["status"] == "ok"
    assert 200 < j["iterations"] < 1700                  # before the fake reaches 1e-12 on its own
    assert max(j["residuals_final"].values()) <= 1e-6
