"""The study design and the driver's task logic, with OpenFOAM replaced by a stand-in."""
import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import driver  # noqa: E402
from cyl import case, study  # noqa: E402


def test_design_sizes_and_tags():
    assert len(study.MEMBERS) == 13 and study.MEMBERS[0] == 10.0 and study.MEMBERS[-1] == 40.0
    assert len(study.GRID) == 840
    grid = study.tasks_grid()
    assert len(grid) == 91 and all(len(t.runs) == 120 for t in grid)
    tags = [r.tag for t in grid for r in t.runs]
    assert len(tags) == len(set(tags)) == 13 * 840
    ref = study.tasks_ref()
    assert len(ref) == 39 and all(r.ref and r.max_iter == 20000 for t in ref for r in t.runs)
    assert not set(tags) & {r.tag for t in ref for r in t.runs}
    # the check runs are grid runs: their results are part of the study
    assert {r.tag for r in study.CHECK_RUNS} <= set(tags)
    assert study.default_start() == dict(n_arc=24, res_ctrl=1e-5, lin_tol=1e-6, n_coarsest=250)
    for s in (study.TUTORIAL, study.CAUTIOUS, study.default_start()):
        assert study.setting_of(study.index_of(s)) == s
    assert study.run_at(12.5, study.default_start()).tag == "Re12.5_n24_rc1e-05_lt1e-06_nc250"


def test_inflow_gives_the_reynolds_number():
    for re_ in study.MEMBERS:
        ubar = 2.0 * case.u_max(re_) / 3.0
        assert abs(ubar * case.D / case.NU - re_) < 1e-12
    assert abs(case.u_max(20.0) - 0.3) < 1e-15            # the article's 2D-1


def expand(arr):
    out = []
    for part in arr.split(","):
        a, _, b = part.partition("-")
        out.extend(range(int(a), int(b or a) + 1))
    return out


class FakeMesh:
    made = []

    def __init__(self, root, re_, na):
        if os.environ.get("FAKE_MESH_FAIL"):
            raise case.OpenFOAMError("blockMesh failed (fake)")
        self.re, self.na, self.cells, self.mesh_ok, self.flow_rate = re_, na, 100 * na, True, 1.0
        FakeMesh.made.append((re_, na))

    def remove(self):
        pass


def fake_run(mesh, run_dir, rc, lt, nc, max_iter=2000):
    conv = rc >= 1e-6
    return dict(status="ok", converged=conv, iterations=10 if conv else max_iter, cpu=0.1 * mesh.na,
                wall=0.1, cells=mesh.cells, cells_iters=mesh.cells * 10, returncode=0,
                residuals_final={}, cd=5.5 + rc, cl=0.01, dp=0.1)


@pytest.fixture
def work(tmp_path, monkeypatch):
    for name, sub in (("RUNS", "runs"), ("TASKS", "tasks"), ("CASES", "cases")):
        monkeypatch.setattr(driver, name, str(tmp_path / sub))
    monkeypatch.setattr(case, "Mesh", FakeMesh)
    monkeypatch.setattr(case, "run", fake_run)
    monkeypatch.setattr(case, "openfoam_version", lambda: "v2606")
    FakeMesh.made = []
    return tmp_path


def test_task_runs_every_setting_once(work):
    t = study.tasks_grid()[4]                      # Re 10, n_arc 32
    assert driver.run_task(t) == 0
    assert FakeMesh.made == [(10.0, 32)]
    files = os.listdir(work / "runs")
    assert len(files) == 120
    x = json.load(open(work / "runs" / f"{t.runs[0].tag}.json"))
    assert x["tag"] == t.runs[0].tag and x["openfoam"] == "v2606" and x["task"] == t.name
    assert x["setting"] == t.runs[0].setting() and x["max_iter"] == 2000
    assert driver.run_task(t) == 0 and FakeMesh.made == [(10.0, 32)]      # nothing left to run
    assert driver.task_done(t) and not os.path.exists(driver._lock(t))


def test_errors_are_recorded_and_given_up_after_three(work, monkeypatch):
    monkeypatch.setenv("FAKE_MESH_FAIL", "1")
    t = study.tasks_ref()[0]
    assert [driver.run_task(t) for _ in range(3)] == [2, 2, 2]
    assert len(driver.errors(t)) == 3 and driver.given_up(t)
    assert driver.run_task(t) == 4
    assert "blockMesh failed" in open(work / "tasks" / driver.errors(t)[0]).read()
    submitted = [i for g in driver.groups(study.tasks_ref()) for i in expand(g.split()[2])]
    assert 0 not in submitted and len(submitted) == 38                   # task 0 not resubmitted


def test_lock_of_a_live_process_blocks_and_a_dead_one_does_not(work):
    t = study.tasks_check()[0]
    os.makedirs(work / "tasks", exist_ok=True)
    driver.write_json(driver._lock(t), dict(host=driver.socket.gethostname(), pid=os.getpid()))
    assert driver.run_task(t) == 5
    driver.write_json(driver._lock(t), dict(host=driver.socket.gethostname(), pid=2 ** 22 + 12345))
    assert driver.run_task(t) == 0


def test_stop_flag_interrupts_between_runs(work, monkeypatch):
    t = study.tasks_grid()[0]
    calls = {"n": 0}

    def run_then_stop(*a, **k):
        calls["n"] += 1
        if calls["n"] == 5:
            driver._STOP["flag"] = True
        return fake_run(*a, **k)

    monkeypatch.setattr(case, "run", run_then_stop)
    driver._STOP["flag"] = False
    try:
        assert driver.run_task(t) == 3
    finally:
        driver._STOP["flag"] = False
    assert len(os.listdir(work / "runs")) == 5
    assert driver.run_task(t) == 0 and len(os.listdir(work / "runs")) == 120


def test_groups_by_resource_class(work):
    g = driver.groups(study.tasks_grid())
    rows = [x.split() for x in g]
    assert {(r[0], r[1]) for r in rows} == {("04:00:00", "2G"), ("12:00:00", "3G")}
    assert sorted(i for r in rows for i in expand(r[2])) == list(range(91))
    assert driver.compress([0, 1, 2, 5, 7, 8]) == "0-2,5,7-8"
