"""The job-file layer: loading, validation, the command runner, caching, the CLI."""
import json
import os
import sys
import textwrap

import pytest

from errbudget import load_job, run_job, JobSpec
from errbudget.jobfile import check_job, init_job
from errbudget.runner import (CommandSolver, JobError, RunFailed, extract_number, render,
                          format_value)
from errbudget.cli import main

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(HERE, "..", "examples", "04_command_line")
PY = sys.executable


# ------------------------------------------------------------------ rendering
def test_placeholders_are_rendered_and_unknown_ones_refused():
    assert render("run --n {{ n }} --dt {{dt}} --x {{x:.2f}}", {"n": 32, "dt": 1e-6, "x": 0.5}) \
        == "run --n 32 --dt 1e-06 --x 0.50"
    with pytest.raises(JobError, match="unknown placeholder"):
        render("{{nope}}", {"n": 1})
    # OpenFOAM-style braces are left alone
    assert render("solvers { p { tolerance {{tol}}; } }", {"tol": 1e-8}) == \
        "solvers { p { tolerance 1e-08; } }"


def test_format_value():
    assert format_value(3) == "3" and format_value(1e-6) == "1e-06" and format_value(0.05) == "0.05"
    assert format_value(True) == "true" and format_value("abc") == "abc"


# ------------------------------------------------------------------ extraction
def test_extract_number_methods():
    assert extract_number({"regex": r"Cd\s*=\s*([-0-9.e+]+)"}, "Cd = 1.5\nCd = 2.5e-1\n") == 0.25
    assert extract_number({"regex": r"Cd\s*=\s*([-0-9.e+]+)", "which": "first"}, "Cd = 1.5\nCd = 0.25") == 1.5
    assert extract_number({"key": "forces.drag"}, '{"forces": {"drag": 3.0}}') == 3.0
    assert extract_number({"key": "rows[1].v"}, '{"rows": [{"v": 1}, {"v": 2}]}') == 2.0
    assert extract_number({"column": 2}, "# t Cd Cl\n1 0.1 0.2\n2 0.3 0.4\n") == 0.3
    assert extract_number({}, "energy -73671.37 kJ/mol\n") == -73671.37
    assert extract_number({"regex": r"E = ([-0-9.D+]+)"}, "E = -1.5D-03") == -1.5e-3
    with pytest.raises(JobError, match="matched nothing"):
        extract_number({"regex": "zzz"}, "nothing")
    with pytest.raises(JobError):
        extract_number({"key": "a.b"}, '{"a": {}}')


# ------------------------------------------------------------------ loading
def _write(tmp_path, text, name="job.yaml"):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text))
    return str(p)


def test_yaml_numbers_without_a_dot_are_numbers(tmp_path):
    pytest.importorskip("yaml")
    p = _write(tmp_path, """
        command: echo {{n}} {{tol}}
        parameters:
          - {name: n, values: [8, 16, 32], start: 16}
          - {name: tol, kind: nested-solver, values: [1e-2, 1e-4, 1e-6]}
        observable: {name: v}
        target: 1e-3
        """)
    job = load_job(p)
    assert job.target == 1e-3
    assert job.parameters[1]["values"] == [1e-2, 1e-4, 1e-6]
    assert job.start == {"n": 16}
    assert job.to_parameters()[1].kind == "nested-solver"


def test_json_and_toml_job_files(tmp_path):
    d = {"command": "echo {{n}}", "parameters": [{"name": "n", "values": [1, 2, 3]}],
         "observable": {"name": "v"}, "target": 0.1}
    p = tmp_path / "job.json"
    p.write_text(json.dumps(d))
    assert load_job(str(p)).target == 0.1
    if sys.version_info >= (3, 11):
        t = tmp_path / "job.toml"
        t.write_text('command = "echo {{n}}"\ntarget = 0.1\n[observable]\nname = "v"\n'
                     '[[parameters]]\nname = "n"\nvalues = [1, 2, 3]\n')
        assert load_job(str(t)).parameters[0]["values"] == [1, 2, 3]


def test_validation_lists_every_problem(tmp_path):
    d = {"command": "run {{n}} {{missing}}", "target": -1,
         "parameters": [{"name": "n", "kind": "weird", "values": [1]},
                        {"name": "n", "values": [1, 2], "start": 7}],
         "observable": "v", "bogus": 1}
    with pytest.raises(JobError) as e:
        JobSpec.from_dict(d)
    msg = str(e.value)
    for frag in ("kind must be", "at least two", "listed twice", "start value 7",
                 "must be positive", "{{missing}}", "unknown field"):
        assert frag in msg, frag


def test_solver_and_command_are_exclusive():
    base = {"parameters": [{"name": "n", "values": [1, 2]}], "observable": "v", "target": 1}
    with pytest.raises(JobError, match="either 'command'"):
        JobSpec.from_dict(dict(base))
    with pytest.raises(JobError, match="not both"):
        JobSpec.from_dict(dict(base, command="x", solver="m:f"))


def test_kind_aliases_and_configurations():
    d = {"command": "x {{n}}", "parameters": [{"name": "n", "kind": "mesh", "values": [1, 2]}],
         "observable": "v", "target": 1, "options": {"configuration": "v1.1"}}
    job = JobSpec.from_dict(d)
    assert job.parameters[0]["kind"] == "discretisation"
    kw = job.tune_kwargs()
    assert kw["joint_check"] is False and kw["cancel_rule"] == "credit"
    with pytest.raises(JobError, match="configuration"):
        JobSpec.from_dict(dict(d, options={"configuration": "v9"}))


# ------------------------------------------------------------------ the runner
SOLVER = """
import sys, json
n = int(sys.argv[1]); tol = float(sys.argv[2])
if n == 999: raise SystemExit("boom")
if tol < 1e-9: print("did not converge")
print(json.dumps({"value": 1.0 + 1.0 / n + tol, "work": n / tol}))
"""


def _job(tmp_path, **extra):
    (tmp_path / "s.py").write_text(SOLVER)
    d = {"command": f"{PY} {{{{job_dir}}}}/s.py {{{{n}}}} {{{{tol}}}}",
         "parameters": [{"name": "n", "values": [2, 4, 8, 16, 32, 999], "start": 8},
                        {"name": "tol", "kind": "nested-solver", "values": [1e-1, 1e-3, 1e-6, 1e-10]}],
         "observable": {"name": "value", "key": "value"}, "cost": {"key": "work"},
         "target": 0.3, "workdir": "runs", "keep_runs": True,
         "failure": {"unconverged_if_regex": "did not converge"}}
    d.update(extra)
    return JobSpec.from_dict(d, base_dir=str(tmp_path))


def test_command_solver_runs_caches_and_records(tmp_path):
    job = _job(tmp_path)
    s = CommandSolver(job, log=None)
    out, cost = s(n=8, tol=1e-3)
    assert out == {"value": pytest.approx(1.126), "converged": True} and cost == pytest.approx(8000)
    assert os.path.exists(os.path.join(job.workdir, "run_0001", "stdout.txt"))
    assert json.load(open(os.path.join(job.workdir, "run_0001", "params.json"))) == {"n": 8, "tol": 1e-3}
    # second call: from the log, no new run directory
    s(n=8, tol=1e-3)
    assert s.n_cached == 1 and not os.path.exists(os.path.join(job.workdir, "run_0002"))
    # a crash is a failed run, recorded as such and replayed as such
    with pytest.raises(RunFailed, match="exit code"):
        s(n=999, tol=1e-3)
    s2 = CommandSolver(job, log=None)            # a new process resumes from runs.jsonl
    assert len(s2.cache) == 2
    with pytest.raises(RunFailed):
        s2(n=999, tol=1e-3)
    assert s2.n_cached == 1
    # not converged: a value is returned but flagged
    out, _ = s2(n=4, tol=1e-10)
    assert out["converged"] is False
    rows = [json.loads(l) for l in open(os.path.join(job.workdir, "runs.jsonl"))]
    assert [r["status"] for r in rows] == ["ok", "failed", "unconverged"]


def test_fresh_ignores_the_log(tmp_path):
    job = _job(tmp_path)
    CommandSolver(job, log=None)(n=8, tol=1e-3)
    s = CommandSolver(job, log=None, fresh=True)
    s(n=8, tol=1e-3)
    assert s.n_cached == 0


def test_template_dir_and_file_observable(tmp_path):
    case = tmp_path / "case"
    (case / "system").mkdir(parents=True)
    (case / "system" / "controlDict").write_text("deltaT {{dt}};\nwriteInterval 10;\n")
    (case / "run.sh").write_text("cat system/controlDict > out.txt; echo 'value {{n}}' >> out.txt\n")
    d = {"command": "sh run.sh", "template_dir": "case",
         "parameters": [{"name": "n", "values": [1, 2, 3]},
                        {"name": "dt", "values": [1e-2, 1e-3]}],
         "observable": {"name": "value", "source": "file", "path": "out.txt",
                        "regex": r"value (\d+)"}, "target": 1}
    job = JobSpec.from_dict(d, base_dir=str(tmp_path))
    s = CommandSolver(job, log=None)
    out, _ = s(n=2, dt=1e-3)
    assert out["value"] == 2.0
    assert open(os.path.join(job.workdir, "run_0001", "system", "controlDict")).read() \
        == "deltaT 0.001;\nwriteInterval 10;\n"


def test_python_solver_in_a_job_file(tmp_path):
    (tmp_path / "mysolver.py").write_text(
        "def run(n, tol):\n    return {'v': 1 + 1/n + tol, 'converged': tol > 1e-9}, n/tol\n")
    d = {"solver": "mysolver.py:run",
         "parameters": [{"name": "n", "values": [2, 4, 8, 16, 32], "start": 8},
                        {"name": "tol", "kind": "nested-solver", "values": [1e-1, 1e-3, 1e-6, 1e-10]}],
         "observable": "v", "target": 0.3, "workdir": "runs"}
    job = JobSpec.from_dict(d, base_dir=str(tmp_path))
    r = run_job(job, log=None, write=False)
    assert r.certified and r.settings["n"] >= 4


# ------------------------------------------------------------------ end to end
def test_run_job_end_to_end_and_resume(tmp_path):
    job = _job(tmp_path, output="res")
    r = run_job(job, log=None)
    assert r.certified and abs(r.error_estimate) <= 0.3
    assert os.path.exists(str(tmp_path / "res.json")) and os.path.exists(str(tmp_path / "res.md"))
    n_runs = r.runs_used
    # running again replays everything from the log: identical answer, no new run
    s = CommandSolver(job, log=None)
    before = len(s.cache)
    r2 = run_job(job, log=None, write=False)
    assert r2.settings == r.settings and len(CommandSolver(job, log=None).cache) == before
    assert n_runs == r2.runs_used


def test_result_round_trip(tmp_path):
    from errbudget import TuneResult
    job = _job(tmp_path, output="res")
    r = run_job(job, log=None)
    back = TuneResult.load(str(tmp_path / "res.json"))
    assert back.settings == r.settings and back.certified == r.certified
    assert str(back) == str(r)
    assert "| `n` |" in back.markdown()


def test_cli_init_check_run_report(tmp_path, capsys, monkeypatch):
    pytest.importorskip("yaml")
    monkeypatch.chdir(tmp_path)
    assert main(["init", "job.yaml"]) == 0
    assert os.path.exists("job.yaml")
    assert main(["init", "job.yaml"]) == 1                   # refuses to overwrite
    # point the template's command at a real solver
    import shutil
    shutil.copy(os.path.join(EXAMPLE, "solver.py"), "solver.py")
    text = open("job.yaml").read().replace("command: python3 {{job_dir}}/solver.py", f"command: {PY} {{{{job_dir}}}}/solver.py")
    open("job.yaml", "w").write(text)
    assert main(["check", "job.yaml", "--no-run"]) == 0
    assert main(["check", "job.yaml"]) == 0
    out = capsys.readouterr().out
    assert "u_final = 0.54030" in out and "the job file works" in out
    assert main(["run", "job.yaml", "--quiet"]) == 0
    out = capsys.readouterr().out
    assert "Target 0.0001 MET" in out and os.path.exists("errbudget_result.json")
    assert main(["report", "errbudget_result.json"]) == 0
    assert main(["report", "errbudget_result.json", "--markdown"]) == 0
    assert main(["runs", "job.yaml", "--last", "3"]) == 0
    assert "steps" in capsys.readouterr().out
    assert main(["version"]) == 0
    assert main([]) == 2


def test_example_job_file_loads():
    pytest.importorskip("yaml")
    job = load_job(os.path.join(EXAMPLE, "job.yaml"))
    assert [p["name"] for p in job.parameters] == ["steps", "tol", "refresh"]
    assert job.cost == {"regex": r"work\s*=\s*([0-9]+)"}


# ------------------------------------------------------------------ after review
def test_check_retries_a_failed_start_after_the_command_is_fixed(tmp_path, capsys):
    d = {"command": "definitely-not-a-program {{n}}", "parameters": [{"name": "n", "values": [1, 2]}],
         "observable": "v", "target": 1, "workdir": "runs"}
    assert check_job(JobSpec.from_dict(d, base_dir=str(tmp_path))) is False
    out = capsys.readouterr().out
    assert "FAILED" in out and "stdout.txt" in out
    d["command"] = "echo {{n}}"
    assert check_job(JobSpec.from_dict(d, base_dir=str(tmp_path))) is True
    assert "v = 2" in capsys.readouterr().out        # the start is the middle value, n=2


def test_changed_job_file_does_not_reuse_old_runs(tmp_path):
    (tmp_path / "s.py").write_text("import sys; print('a =', sys.argv[1]); print('b =', 2*int(sys.argv[1]))")
    d = {"command": f"{PY} {{{{job_dir}}}}/s.py {{{{n}}}}", "parameters": [{"name": "n", "values": [1, 2, 3]}],
         "observable": {"name": "v", "regex": r"a = (\d+)"}, "target": 1, "workdir": "runs"}
    job = JobSpec.from_dict(d, base_dir=str(tmp_path))
    assert CommandSolver(job, log=None)(n=1)[0]["v"] == 1.0
    d["observable"] = {"name": "v", "regex": r"b = (\d+)"}
    job2 = JobSpec.from_dict(d, base_dir=str(tmp_path))
    s = CommandSolver(job2, log=None)
    assert s.n_stale == 1 and s(n=1)[0]["v"] == 2.0 and s.n_cached == 0
    # run ids keep counting across versions and with --fresh
    s3 = CommandSolver(job2, log=None, fresh=True)
    s3(n=1)
    ids = [r["run_id"] for r in __import__("errbudget.runner", fromlist=["_read_log"])._read_log(s.log_path)]
    assert ids == ["0001", "0002", "0003"]


def test_failure_patterns_see_redirected_output(tmp_path):
    (tmp_path / "s.py").write_text("import sys; print('value 1\\nSIMPLE solution converged' if sys.argv[1]=='1' else 'value 2')")
    d = {"command": f"{PY} {{{{job_dir}}}}/s.py {{{{n}}}} > log.txt", "parameters": [{"name": "n", "values": [1, 2, 3]}],
         "observable": {"name": "v", "source": "file", "path": "log.txt", "regex": r"value (\d+)"},
         "target": 1, "workdir": "runs", "failure": {"converged_if_regex": "SIMPLE solution converged"}}
    s = CommandSolver(JobSpec.from_dict(d, base_dir=str(tmp_path)), log=None)
    assert s(n=1)[0]["converged"] is True
    assert s(n=2)[0]["converged"] is False
    # a file that is neither the observable's nor the cost's: failure.path
    d["observable"] = {"name": "v", "regex": r"value (\d+)"}
    d["command"] = f"{PY} {{{{job_dir}}}}/s.py {{{{n}}}} | tee log.txt"
    d["failure"] = {"failed_if_regex": "SIMPLE", "path": "log.txt"}
    s = CommandSolver(JobSpec.from_dict(d, base_dir=str(tmp_path)), log=None)
    with pytest.raises(RunFailed, match="matched"):
        s(n=1)


def test_timeout_kills_the_whole_process_group(tmp_path):
    (tmp_path / "slow.py").write_text("import time, sys; time.sleep(3); open('done.txt','w').write('x')\n")
    d = {"command": f"echo start && {PY} {{{{job_dir}}}}/slow.py {{{{n}}}}", "timeout": 0.5,
         "parameters": [{"name": "n", "values": [1, 2]}], "observable": "v", "target": 1, "workdir": "runs"}
    job = JobSpec.from_dict(d, base_dir=str(tmp_path))
    with pytest.raises(RunFailed, match="timed out"):
        CommandSolver(job, log=None)(n=1)
    import time
    time.sleep(3.5)
    assert not os.path.exists(os.path.join(job.workdir, "run_0001", "done.txt"))


def test_validation_catches_engine_refusals_and_bad_fields(tmp_path):
    base = {"command": "echo {{r}}", "observable": "v", "target": 1}
    with pytest.raises(JobError, match="refinement direction"):
        JobSpec.from_dict(dict(base, parameters=[{"name": "r", "values": [0.9, 1.0, 1.2]}]))
    job = JobSpec.from_dict(dict(base, parameters=[{"name": "r", "values": [0.9, 1.0, 1.2], "h": "inverse"}]))
    assert job.to_parameters()[0].h[0] > job.to_parameters()[0].h[-1]
    for bad, frag in [
        (dict(base, parameters=[{"name": "r", "values": [1, 2], "h": [1.0]}]), "h must be"),
        (dict(base, parameters=[{"name": "r", "values": [1, 2]}], timeout="1h"), "timeout"),
        (dict(base, parameters=[{"name": "r", "values": [1, 2]}], cost={"source": "file", "regex": "x"}), "needs 'path'"),
        (dict(base, parameters=[{"name": "r", "values": [1, 2]}], observable={"name": "v", "regex": "("}), "does not compile"),
        (dict(base, parameters=[{"name": "r", "values": [1, 2]}], observable={"name": "v", "regex": "x"}), "capture group"),
        (dict(base, parameters=[{"name": "r", "values": [1, 2]}], options={"resid_margn": 1}), "not an option"),
        (dict(base, parameters=[{"name": "r", "values": [1, 2]}], failure={"fail_regex": "x"}), "not a failure rule"),
        ({"solver": "nope.py:run", "parameters": [{"name": "r", "values": [1, 2]}], "observable": "v", "target": 1}, "does not exist"),
    ]:
        with pytest.raises(JobError, match=frag):
            JobSpec.from_dict(bad, base_dir=str(tmp_path))
    job = JobSpec.from_dict(dict(base, parameters=[{"name": "r", "values": [1, 2]}], options={"resid_margin": "1e-3"}))
    assert job.tune_kwargs()["resid_margin"] == 1e-3


def test_null_json_and_fortran_columns(tmp_path):
    with pytest.raises(JobError, match="could not read"):
        extract_number({"key": "flux"}, '{"flux": null}')
    assert extract_number({"column": 2}, "1 -1.5D-03\n") == -1.5e-3
    assert extract_number({"column": 2}, "0.5 (0.12 0.3 0.0)\n") == 0.12
    (tmp_path / "m.py").write_text("def run(n):\n    return None\n")
    job = JobSpec.from_dict({"solver": "m.py:run", "parameters": [{"name": "n", "values": [1, 2]}],
                             "observable": "v", "target": 1, "workdir": "runs"}, base_dir=str(tmp_path))
    with pytest.raises(RunFailed, match="not a number"):
        CommandSolver(job, log=None)(n=1)


def test_output_directory_is_created_and_log_table_survives_a_truncated_line(tmp_path):
    job = _job(tmp_path, output="results/deep/out")
    r = run_job(job, log=None)
    assert os.path.exists(str(tmp_path / "results" / "deep" / "out.json")) and r.certified
    with open(os.path.join(job.workdir, "runs.jsonl"), "a") as f:
        f.write('{"run_id": "9999", "params": {"n": 1')          # a line cut off mid-write
    from errbudget.report import run_table
    assert "9999" not in run_table(os.path.join(job.workdir, "runs.jsonl"))
    s = CommandSolver(job, log=None)
    assert len(s.cache) == r.runs_used


def test_check_renders_templates_before_running(tmp_path, capsys):
    case = tmp_path / "case"
    case.mkdir()
    (case / "input.txt").write_text("n = {{n}}; bad = {{nope}}\n")
    d = {"command": "cat input.txt", "template_dir": "case",
         "parameters": [{"name": "n", "values": [1, 2]}], "observable": "v", "target": 1}
    with pytest.raises(JobError, match="unknown placeholder"):
        check_job(JobSpec.from_dict(d, base_dir=str(tmp_path)), run=False)
