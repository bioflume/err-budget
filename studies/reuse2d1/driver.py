#!/usr/bin/env python3
"""Run and list the tasks of the reuse study (PROTOCOL.md, RUNBOOK.md).

    python3 driver.py --sweep grid --count | --list [--json] | --status | --groups
    python3 driver.py --sweep grid --index K          # run task K (one Slurm array element)
    python3 driver.py --check-report                  # after the check job: CHECK PASSED/FAILED

A task is one mesh at one member: it builds the mesh once and runs every setting of it that
has no stored result, writing one JSON per run to $REUSE_WORK/runs.  A failed or
unconverged run is a result.  An infrastructure error (blockMesh, checkMesh, a crash of the
driver) is written to $REUSE_WORK/tasks/<task>.error-<k>.txt; the task is tried again at the
next submission, at most three times.

Exit codes: 0 done, 2 error recorded, 3 interrupted (time limit; resubmit), 4 given up,
5 running in another job.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cyl import case, study  # noqa: E402

WORK = os.environ.get("REUSE_WORK", os.path.join(HERE, "work"))
RUNS, TASKS, CASES = (os.path.join(WORK, d) for d in ("runs", "tasks", "cases"))
MAX_ATTEMPTS = 3
_STOP = {"flag": False}


def code_version():
    p = os.path.join(HERE, "VERSION")
    return open(p).readline().strip() if os.path.exists(p) else "working tree"


def result_path(run):
    return os.path.join(RUNS, f"{run.tag}.json")


def stored(run):
    return os.path.exists(result_path(run))


def load(run):
    with open(result_path(run)) as f:
        return json.load(f)


def write_json(path, obj):
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
    os.replace(tmp, path)


# ------------------------------------------------------------------ locks and attempts
def _lock(task):
    return os.path.join(TASKS, f"{task.name}.running")


def _alive(info):
    """Is the job or process that holds a lock still running?"""
    job = info.get("slurm_job")
    if job:
        try:
            r = subprocess.run(["squeue", "-h", "-j", str(job)], capture_output=True, text=True, timeout=60)
            return r.returncode == 0 and bool(r.stdout.strip())
        except (OSError, subprocess.SubprocessError):
            return True                       # cannot tell: assume it runs
    if info.get("host") == socket.gethostname():
        try:
            os.kill(int(info.get("pid", -1)), 0)
            return True
        except (OSError, ValueError):
            return False
    return True


def take_lock(task):
    p = _lock(task)
    if os.path.exists(p):
        try:
            info = json.load(open(p))
        except (OSError, ValueError):
            info = {}
        if _alive(info):
            return False
    write_json(p, dict(host=socket.gethostname(), pid=os.getpid(),
                       slurm_job=os.environ.get("SLURM_JOB_ID"), started=time.strftime("%Y-%m-%dT%H:%M:%S")))
    return True


def drop_lock(task):
    try:
        os.remove(_lock(task))
    except OSError:
        pass


def errors(task):
    return sorted(f for f in os.listdir(TASKS) if f.startswith(f"{task.name}.error-")) if os.path.isdir(TASKS) else []


def task_done(task):
    return all(stored(r) for r in task.runs)


def given_up(task):
    return not task_done(task) and len(errors(task)) >= MAX_ATTEMPTS


# ------------------------------------------------------------------ running
def run_task(task, dry_run=False):
    for d in (RUNS, TASKS, CASES):
        os.makedirs(d, exist_ok=True)
    todo = [r for r in task.runs if not stored(r)]
    if not todo:
        print(f"{task.name}: every run has a result")
        return 0
    if given_up(task):
        print(f"{task.name}: given up after {MAX_ATTEMPTS} errors: {', '.join(errors(task))}")
        return 4
    if dry_run:
        print(f"{task.name}: would run {len(todo)} of {len(task.runs)} runs, e.g. {todo[0].tag}")
        return 0
    if not take_lock(task):
        print(f"{task.name}: RUNNING ELSEWHERE (lock {_lock(task)})")
        return 5
    signal.signal(signal.SIGUSR1, lambda *a: _STOP.update(flag=True))
    root = os.path.join(CASES, task.name)
    meta = dict(host=socket.gethostname(), slurm_job=os.environ.get("SLURM_JOB_ID"),
                slurm_array=os.environ.get("SLURM_ARRAY_JOB_ID"), code=code_version())
    try:
        os.makedirs(root, exist_ok=True)
        meta["openfoam"] = case.openfoam_version()
        t0 = time.time()
        mesh = case.Mesh(root, task.re, task.n_arc)
        meta.update(mesh_cells=mesh.cells, mesh_ok=mesh.mesh_ok, mesh_s=round(time.time() - t0, 1),
                    flow_rate=mesh.flow_rate)
        print(f"{task.name}: mesh {mesh.cells} cells (checkMesh {'OK' if mesh.mesh_ok else 'NOT OK'}), "
              f"{len(todo)} runs to do", flush=True)
        for k, r in enumerate(todo):
            if _STOP["flag"]:
                print(f"{task.name}: stopped before {r.tag} (time limit); resubmit")
                return 3
            started = time.strftime("%Y-%m-%dT%H:%M:%S")
            res = case.run(mesh, os.path.join(root, "run"), r.res_ctrl, r.lin_tol, r.n_coarsest,
                           max_iter=r.max_iter)
            out = dict(r.to_dict(), setting=r.setting(), task=task.name, started=started,
                       finished=time.strftime("%Y-%m-%dT%H:%M:%S"), u_max=case.u_max(r.re),
                       max_iter=r.max_iter, **meta, **res)
            write_json(result_path(r), out)
            print(f"  [{k + 1}/{len(todo)}] {r.tag}: {res['status']}, "
                  f"{'converged' if res['converged'] else 'NOT converged'} in {res['iterations']} it, "
                  f"cd {res['cd']:.8g}, cpu {res['cpu']:.1f} s", flush=True)
        mesh.remove()
        shutil.rmtree(root, ignore_errors=True)
        return 0
    except Exception:                                     # noqa: BLE001 - recorded and retried
        k = len(errors(task)) + 1
        p = os.path.join(TASKS, f"{task.name}.error-{k}.txt")
        with open(p, "w") as f:
            f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {json.dumps(meta)}\n{traceback.format_exc()}")
        print(f"{task.name}: ERROR {k} of {MAX_ATTEMPTS} recorded in {p}\n{traceback.format_exc()}")
        return 2
    finally:
        drop_lock(task)


# ------------------------------------------------------------------ listing
def compress(idx):
    """[0, 1, 2, 5, 7, 8] -> '0-2,5,7-8' (a Slurm --array list)."""
    out, idx, i = [], sorted(idx), 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
            j += 1
        out.append(f"{idx[i]}" if i == j else f"{idx[i]}-{idx[j]}")
        i = j + 1
    return ",".join(out)


def groups(tasks):
    """'time mem array' per resource class, for the tasks still to run."""
    by = collections.defaultdict(list)
    for i, t in enumerate(tasks):
        if task_done(t) or given_up(t) or (os.path.exists(_lock(t)) and _alive(json.load(open(_lock(t))))):
            continue
        by[study.resource_class(t)].append(i)
    return [f"{tm} {mem} {compress(ix)}" for (tm, mem), ix in sorted(by.items())]


def status(name, tasks):
    n_runs = sum(len(t.runs) for t in tasks)
    res = [load(r) for t in tasks for r in t.runs if stored(r)]
    st = collections.Counter("converged" if x["converged"] else ("failed" if x["status"] == "failed" else "at the cap")
                             for x in res)
    core_h = sum(float(x.get("cpu") or 0.0) for x in res) / 3600
    done = sum(task_done(t) for t in tasks)
    running = sum(os.path.exists(_lock(t)) for t in tasks)
    gave = [t.name for t in tasks if given_up(t)]
    errs = sum(len(errors(t)) for t in tasks)
    print(f"{name}: tasks {done}/{len(tasks)} done, {running} running, {len(gave)} given up, {errs} error files; "
          f"runs {len(res)}/{n_runs} stored ({', '.join(f'{k} {v}' for k, v in sorted(st.items()))}); "
          f"simpleFoam CPU {core_h:.2f} core-hours")
    for g in gave:
        print(f"  GIVEN UP: {g}")


def check_report():
    """PROTOCOL.md 11: the three check runs converge and c_D at Re 20, n_arc 12 is within
    3 % of the benchmark; the OpenFOAM release is v2606."""
    ok = True
    for r in study.CHECK_RUNS:
        if not stored(r):
            print(f"{r.tag}: no result")
            ok = False
            continue
        x = load(r)
        print(f"{r.tag}: {x['status']}, {'converged' if x['converged'] else 'NOT converged'} in "
              f"{x['iterations']} iterations, cd {x['cd']:.6f}, cpu {x['cpu']:.2f} s, cells {x['cells']}, "
              f"OpenFOAM {x.get('openfoam')}")
        ok &= bool(x["converged"])
        if str(x.get("openfoam")) != "v2606":
            print(f"  version {x.get('openfoam')}, not v2606")
            ok = False
    r = study.CHECK_RUNS[1]
    if stored(r):
        dev = abs(load(r)["cd"] - study.BENCHMARK[20.0]) / study.BENCHMARK[20.0]
        print(f"Re 20, n_arc 12: cd {load(r)['cd']:.6f}, {dev:.2%} from the benchmark {study.BENCHMARK[20.0]}")
        ok &= dev <= 0.03
    for t in study.tasks_check():
        for e in errors(t):
            print(f"error file {e}:\n" + open(os.path.join(TASKS, e)).read()[-1500:])
    print("CHECK PASSED" if ok else "CHECK FAILED")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", choices=sorted(study.SWEEPS))
    ap.add_argument("--index", type=int)
    ap.add_argument("--count", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--groups", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check-report", action="store_true")
    a = ap.parse_args(argv)
    if a.check_report:
        return check_report()
    if not a.sweep:
        ap.error("--sweep is needed")
    tasks = study.sweep(a.sweep)
    if a.count:
        print(len(tasks))
    elif a.list:
        for i, t in enumerate(tasks):
            row = dict(index=i, task=t.name, re=t.re, n_arc=t.n_arc, runs=len(t.runs),
                       stored=sum(stored(r) for r in t.runs), time_mem=study.resource_class(t))
            print(json.dumps(row) if a.json else f"{i:4d} {t.name:24s} {row['stored']}/{row['runs']} stored")
    elif a.status:
        status(a.sweep, tasks)
    elif a.groups:
        print("\n".join(groups(tasks)))
    elif a.index is not None:
        return run_task(tasks[a.index], dry_run=a.dry_run)
    else:
        ap.error("nothing to do")
    return 0


if __name__ == "__main__":
    sys.exit(main())
