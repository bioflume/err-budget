#!/usr/bin/env python3
"""Run driver of the airfoil study (PROTOCOL.md, RUNBOOK.md).

    python3 driver.py --sweep stage1 --list            # index, tag, cores, cells, status
    python3 driver.py --sweep stage1 --count           # number of runs (= array size)
    python3 driver.py --sweep stage1 --index 17        # run index 17 (skips it if done)
    python3 driver.py --sweep stage1 --index 17 --dry-run [--case-dir DIR]
                                                       # write the OpenFOAM case only
    python3 driver.py --sweep stage1 --status          # done / failed / errors / missing
    python3 driver.py --sweep stage1 --groups          # unfinished indices by resource class
    python3 driver.py --run-setting m=0.06 t=0.08 n_surf=128 dy1=2e-5 radius=500 \\
                      res_ctrl=1e-5 gamg_tol=1e-7 --print
                                                       # one run, for the errbudget job file

Every run has a deterministic tag (study.Run.tag).  Results go to one store,
$AIRFOIL_WORK/runs/<tag>.json, shared by every sweep: a run that one stage has done
is not done again by another.  A run is a result when OpenFOAM ran it to an end:
converged, stopped at the iteration cap, stopped as stalled, or diverged after it
started (status "failed").  Anything else (a missing module, a mesh that checkMesh
rejects, a crash before the first iteration) is an error: it is written to
<tag>.error-<k>.json and the run is tried again at the next submission, at most
three times.  A run interrupted by the Slurm time limit (SIGUSR1) stops cleanly with
its fields written and resumes from there at the next submission.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import signal
import socket
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402

from af import foamcase, foamio, study  # noqa: E402
from af.cgrid import MESH_VERSION, MeshSpec, build  # noqa: E402
from af.naca import Aerofoil  # noqa: E402

WORK = os.environ.get("AIRFOIL_WORK") or os.path.join("/nobackup", os.environ.get("USER", "user"), "airfoil")
RUNS = os.path.join(WORK, "runs")
CASES = os.path.join(WORK, "cases")
MAX_ERRORS = 3
POLL_S = float(os.environ.get("AIRFOIL_POLL_S", 10.0))
MPIRUN = os.environ.get("AIRFOIL_MPIRUN", "mpirun")
MPIRUN_OPTS = os.environ.get("AIRFOIL_MPIRUN_OPTS", "--bind-to none --oversubscribe").split()

_STOP = {"flag": False}


def _on_usr1(signum, frame):
    _STOP["flag"] = True


def code_hash() -> str:
    h = hashlib.sha1()
    for rel in ("driver.py", "af/naca.py", "af/grading.py", "af/cgrid.py", "af/foamcase.py",
                "af/foamio.py", "af/study.py"):
        with open(os.path.join(HERE, rel), "rb") as fh:
            h.update(fh.read())
    return h.hexdigest()[:12]


def now() -> str:
    return dt.datetime.now().isoformat(timespec="seconds")


# ---------------------------------------------------------------- store
def result_path(tag):
    return os.path.join(RUNS, f"{tag}.json")


def error_paths(tag):
    if not os.path.isdir(RUNS):
        return []
    return sorted(os.path.join(RUNS, f) for f in os.listdir(RUNS)
                  if f.startswith(tag + ".error-") and f.endswith(".json"))


def is_stale(path) -> bool:
    """A result made with an older mesh (no mesh_version: the first generator) or an older
    case set-up (no setup_version: as registered, before A3)."""
    with open(path) as fh:
        r = json.load(fh)
    return r.get("mesh_version", 1) != MESH_VERSION or r.get("setup_version", 1) != study.SETUP_VERSION


def set_aside_stale(tag) -> bool:
    """Move a stale result (older mesh or set-up) and its log to runs/stale/ (kept, never deleted)."""
    p = result_path(tag)
    if not (os.path.exists(p) and is_stale(p)):
        return False
    d = os.path.join(RUNS, "stale")
    os.makedirs(d, exist_ok=True)
    for ext in (".json", ".log.txt"):
        src = os.path.join(RUNS, tag + ext)
        if os.path.exists(src):
            os.replace(src, os.path.join(d, tag + ext))
    print(f"set aside a result made with an older mesh or set-up: {tag} -> {d}/", flush=True)
    return True


def status_of(run) -> str:
    p = result_path(run.tag)
    if os.path.exists(p):
        if is_stale(p):
            return "stale"
        return "failed" if json.load(open(p)).get("status") == "failed" else "done"
    if len(error_paths(run.tag)) >= MAX_ERRORS:
        return "gave-up"
    if running_elsewhere(run):
        return "running"
    if error_paths(run.tag):
        return "error"
    if os.path.exists(os.path.join(CASES, run.tag, "STATE.json")):
        return "started"
    return "missing"


def _atomic_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
    os.replace(tmp, path)


def write_error(run, stage, message, case=None):
    k = len(error_paths(run.tag)) + 1
    tail = ""
    if case:
        for name in ("log.checkMesh", "log.simpleFoam.0", "log.decomposePar"):
            p = os.path.join(case, name)
            if os.path.exists(p):
                tail += f"\n==> {name} (tail)\n" + "".join(open(p, errors="replace").readlines()[-40:])
    _atomic_json(os.path.join(RUNS, f"{run.tag}.error-{k}.json"),
                 dict(tag=run.tag, stage=stage, message=message, log_tail=tail[-8000:],
                      host=socket.gethostname(), time=now(), slurm=slurm_ids(), code=code_hash()))
    print(f"ERROR ({stage}) {run.tag}: {message}", flush=True)


def slurm_ids():
    return {k: os.environ.get(v) for k, v in (("job", "SLURM_JOB_ID"), ("array_job", "SLURM_ARRAY_JOB_ID"),
                                              ("task", "SLURM_ARRAY_TASK_ID"), ("cpus", "SLURM_CPUS_PER_TASK"))}


# ---------------------------------------------------------------- one run at a time
# A run is executed by one process at a time.  Two copies of the same run (a sweep
# submitted again while its first array is still running, 10 Nov 2026) share the case
# directory: their logs interleave and the first to finish deletes the case under the
# other.  Each execution therefore holds runs/<tag>.running, created exclusively; a lock
# left by a job that no longer exists (killed at its time limit) is taken over.
def lock_path(tag):
    return os.path.join(RUNS, f"{tag}.running")


def _slurm_name(ids):
    if ids.get("array_job") and ids.get("task") is not None:
        return f"{ids['array_job']}_{ids['task']}"
    return ids.get("job")


def lock_alive(info) -> bool:
    """Whether the process that wrote this lock may still be running."""
    mine = _slurm_name(slurm_ids())
    theirs = _slurm_name(info.get("slurm") or {})
    if theirs:
        if theirs == mine:
            return False                        # this very job, requeued: its lock is stale
        if shutil.which("squeue"):
            r = subprocess.run(["squeue", "-h", "-j", theirs, "-o", "%T"], capture_output=True, text=True)
            if r.returncode == 0:
                return bool(r.stdout.strip())
            return "Invalid job id" not in (r.stderr or "")
    if info.get("host") == socket.gethostname():
        try:
            os.kill(int(info.get("pid", -1)), 0)
            return True
        except (OSError, ValueError):
            return False
    return True                                 # elsewhere and unverifiable: assume it runs


def read_lock(tag):
    try:
        with open(lock_path(tag)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def running_elsewhere(run) -> bool:
    info = read_lock(run.tag)
    return bool(info) and lock_alive(info)


def acquire_lock(tag) -> bool:
    os.makedirs(RUNS, exist_ok=True)
    me = dict(host=socket.gethostname(), pid=os.getpid(), slurm=slurm_ids(), time=now())
    for _ in range(2):
        try:
            fd = os.open(lock_path(tag), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            info = read_lock(tag)
            if info is not None and lock_alive(info):
                return False
            if read_lock(tag) == info:          # unchanged since judged stale: take it over
                try:
                    os.remove(lock_path(tag))
                except FileNotFoundError:
                    pass
            continue
        with os.fdopen(fd, "w") as fh:
            json.dump(me, fh)
        return True
    return False


def release_lock(tag):
    info = read_lock(tag)
    if info and info.get("pid") == os.getpid() and info.get("host") == socket.gethostname():
        try:
            os.remove(lock_path(tag))
        except FileNotFoundError:
            pass


# ---------------------------------------------------------------- one run
def make_grid(run):
    foil = Aerofoil(run.m, run.t, p=study.P_CAMBER, alpha_deg=study.ALPHA_DEG)
    return build(foil, MeshSpec(run.n_surf, run.dy1, run.radius))


def write_case(run, case, diag=False):
    """The case and its mesh, after the generator's own check."""
    grid = make_grid(run)
    chk = grid.check()
    if chk["inverted"] or chk["faces_over_70_deg"]:
        raise RuntimeError(f"generator check: {chk['inverted']} inverted cells, "
                           f"{chk['faces_over_70_deg']} faces above 70 degrees")
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)
    foamcase.write_case(case, run, grid, run.num, diag=diag)
    return grid, chk


DIAG_SWEEPS = ("numerics",)       # sweeps whose runs keep two late times and their fields
# the settling measurement: no stall watchdog (only the iteration cap stops a run), and the
# whole force history kept in runs/history/<tag>.csv
LONG_SWEEPS = ("settle",)


def field_diagnostics(run, case):
    """For the numerics sweep: from the two latest written times, the flow at the wall
    (first-cell velocity along the surface, pressure coefficient: where the boundary layer
    separates) and where the solution still changes between them.  Saves the fields to
    runs/fields/<tag>.npz for the laptop and returns a summary for the result JSON."""
    grid = make_grid(run)
    nn = grid.n_in + grid.n_out
    n = grid.n_cells()
    times = sorted((float(t), t) for t in os.listdir(case)
                   if foamio.NUM_DIR.fullmatch(t) and float(t) > 0 and os.path.isdir(os.path.join(case, t)))
    if not times:
        raise RuntimeError("no written time")
    (t1, d1), prev = times[-1], (times[-2] if len(times) > 1 else None)

    def field(d, name):
        path = os.path.join(case, d, name)
        return foamio.internal_field(path, n) if os.path.exists(path) else None
    U = field(d1, "U")[:, :2]
    p = field(d1, "p")
    nut = field(d1, "nut")
    turb = {k: field(d1, k) for k in ("k", "omega", "nuTilda")}
    G = grid.G
    a, b, c, d = G[:-1, :-1], G[1:, :-1], G[1:, 1:], G[:-1, 1:]
    area = 0.5 * np.abs((c[..., 0] - a[..., 0]) * (d[..., 1] - b[..., 1])
                        - (d[..., 0] - b[..., 0]) * (c[..., 1] - a[..., 1])).reshape(-1)
    cen = grid.centroids().reshape(-1, 2)
    # the wall: columns n_lw .. n_lw + n_wall - 1, lower side trailing edge -> leading edge,
    # then upper side leading edge -> trailing edge; first cells j = 0
    iw = np.arange(grid.n_lw, grid.n_lw + grid.n_wall)
    cw = iw * nn
    e = G[iw + 1, 0] - G[iw, 0]
    tvec = e / np.linalg.norm(e, axis=1)[:, None]
    upper = iw >= grid.n_lw + grid.n_wall // 2
    u_down = np.einsum("ij,ij->i", U[cw], tvec) * np.where(upper, 1.0, -1.0)
    xw = 0.5 * (G[iw + 1, 0, 0] + G[iw, 0, 0])
    lw = np.linalg.norm(e, axis=1)
    wall = {}
    for side, sel in (("upper", upper), ("lower", ~upper)):
        o = np.argsort(xw[sel])
        x, u, ln = xw[sel][o], u_down[sel][o], lw[sel][o]
        aft = x > 0.05                          # away from the stagnation point
        rev = aft & (u < 0)
        # separation_x: the first wall face (from the leading edge) with reversed flow;
        # reversed_fraction: its share of the surface length beyond x = 0.05
        wall[side] = dict(separation_x=float(x[rev][0]) if rev.any() else None,
                          reversed_fraction=float(ln[rev].sum() / max(ln[aft].sum(), 1e-300)),
                          x=x.tolist(), u_first_cell=u.tolist(), cp=(2 * p[cw][sel][o]).tolist())
    out = dict(time=t1, wall=wall)
    if turb["k"] is not None:
        i = int(np.argmax(turb["k"]))
        out["k_max"] = dict(value=float(turb["k"][i]), x=float(cen[i, 0]), y=float(cen[i, 1]))
    i = int(np.argmax(nut))
    out["nut_max"] = dict(nu_ratio=float(nut[i] / study.NU), x=float(cen[i, 0]), y=float(cen[i, 1]))
    dU = None
    if prev:
        t0, d0 = prev
        U0 = field(d0, "U")[:, :2]
        dU = np.linalg.norm(U - U0, axis=1) / max(t1 - t0, 1.0)
        w = dU ** 2 * area
        col = np.repeat(np.arange(grid.n_a), nn)
        row = np.tile(np.arange(nn), grid.n_a)
        inner = row < grid.n_in
        on_wall = (col >= grid.n_lw) & (col < grid.n_lw + grid.n_wall)
        up_cols = col >= grid.n_lw + grid.n_wall // 2
        regions = {
            "upper surface, x < 0.5": inner & on_wall & up_cols & (cen[:, 0] < 0.5),
            "upper surface, x >= 0.5": inner & on_wall & up_cols & (cen[:, 0] >= 0.5),
            "lower surface": inner & on_wall & ~up_cols,
            "wake (inner layer)": inner & ~on_wall,
            "outer layer": ~inner,
        }
        tot = float(w.sum()) or 1.0
        top = np.argsort(-dU)[:10]
        out["change"] = dict(between=[t0, t1], max_per_iter=float(dU.max()),
                             share=({k: float(w[m].sum() / tot) for k, m in regions.items()}),
                             top=[dict(x=float(cen[i, 0]), y=float(cen[i, 1]), dU=float(dU[i])) for i in top])
    os.makedirs(os.path.join(RUNS, "fields"), exist_ok=True)
    arrays = dict(G=G.astype(np.float32), U=U.astype(np.float32), p=p.astype(np.float32),
                  nut=nut.astype(np.float32), time=np.array([t1]))
    for k, v in turb.items():
        if v is not None:
            arrays[k] = v.astype(np.float32)
    if dU is not None:
        arrays["dU"] = dU.astype(np.float32)
    np.savez_compressed(os.path.join(RUNS, "fields", f"{run.tag}.npz"), **arrays)
    return out


def check_mesh(case, grid):
    """checkMesh on the written mesh; a failed check other than the aspect ratio is an error."""
    text = sh(["checkMesh"], case, "log.checkMesh", check=False)
    cm = foamio.parse_check_mesh(text)
    fails = foamio.check_mesh_failures(text)
    if not cm.get("cells"):
        fails.append("checkMesh did not report a cell count")
    elif cm["cells"] != grid.n_cells():
        fails.append(f"checkMesh counts {cm['cells']} cells, the generator {grid.n_cells()}")
    cm["failures"] = fails
    return cm


def sh(cmd, case, log, check=True):
    with open(os.path.join(case, log), "w") as fo:
        r = subprocess.run(cmd, cwd=case, stdout=fo, stderr=subprocess.STDOUT)
    if check and r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} exited with {r.returncode} (see {log})")
    return open(os.path.join(case, log), errors="replace").read()


def solver_cmd(run):
    if run.np > 1:
        return [MPIRUN, "-np", str(run.np)] + MPIRUN_OPTS + ["simpleFoam", "-parallel"]
    return ["simpleFoam"]


class SequenceInterrupted(Exception):
    """The job was asked to stop during the first run of a sequenced start."""


def sequence_start(run, case):
    """A4 (PROTOCOL.md 5.1, 14): a run finer than the Stage-1 lists starts from the solution
    of the same design on the coarse SEQUENCE mesh of its own radius.  That first run is
    solved in case/seq (from the free stream, with the watchdog), mapped onto the run's mesh
    with mapFields, and its cost recorded; a first run that diverges is an error."""
    src = study.sequence_run(run)
    sdir = os.path.join(case, "seq")
    grid_s, _ = write_case(src, sdir)
    cm_s = check_mesh(sdir, grid_s)
    if cm_s["failures"]:
        raise RuntimeError(f"checkMesh of the first run {src.tag}: " + "; ".join(cm_s["failures"]))
    qc = grid_s.foil.quarter_chord()
    tmp = dict(segments=0, resume_from=0, diag=False)
    rc, st, reason = run_segment(src, sdir, tmp, [float(qc[0]), float(qc[1])])
    log = os.path.join(sdir, "log.simpleFoam.0")
    if reason == "interrupted":
        raise SequenceInterrupted()
    kind = foamio.log_failure_kind(log) if rc != 0 else None
    if reason == "diverged" or st.nonfinite or kind == "numerical":
        raise RuntimeError(f"the first run {src.tag} diverged after {st.it} iterations")
    if rc != 0 and not (reason in ("stalled", "settled") and kind == "io"):
        raise RuntimeError(f"the first run {src.tag}: simpleFoam exited with {rc} ({kind or 'no message'})")
    latest = foamio.latest_time(sdir)
    if not latest or float(latest) <= 0:
        raise RuntimeError(f"the first run {src.tag} wrote no solution")
    foamcase.write_map_fields_dict(case)
    sh(["mapFields", "seq", "-sourceTime", "latestTime"], case, "log.mapFields")
    shutil.copy(log, os.path.join(case, "log.sequence.simpleFoam"))
    # its final L/D: with the far-field run at R = 500 it gives the correction of A5
    try:
        _, c = foamio.read_coefficients(sdir)
        ld_first = float(c["Cl"][-1] / c["Cd"][-1])
    except Exception:
        ld_first = math.nan
    settled = bool(st.converged) or reason == "settled"
    out = dict(tag=src.tag, iterations=int(st.it), converged=settled, ld=ld_first,
               stop="settled" if reason == "settled" else ("converged" if st.converged else (reason or "endTime")),
               cells=int(grid_s.n_cells()), cells_iters=float(grid_s.n_cells()) * st.it,
               cpu_s=float(st.exec_s), wall_s=float(tmp.get("wall_s", 0.0)),
               residuals_final=st.final_res, time=float(latest))
    shutil.rmtree(sdir, ignore_errors=True)
    print(f"first run {src.tag}: {out['stop']} after {out['iterations']} iterations, "
          f"mapped onto {run.tag}", flush=True)
    return out


def run_segment(run, case, state, cofr):
    """One solver segment with the watchdog.  Returns (returncode, LogState, reason)."""
    seg = state["segments"]
    log = os.path.join(case, f"log.simpleFoam.{seg}")
    start = state.get("resume_from", 0)
    diag = state.get("diag", False)
    foamcase.set_control(case, run, stop="endTime", cofr=cofr, diag=diag)
    st = foamio.LogState(10**9 if state.get("long") else study.STALL_ITER, first_iter=start)
    tail = foamio.CoeffTail(case, study.SETTLE_WINDOW) if run.settle else None
    reason = None
    t0 = time.time()
    with open(log, "w") as fo:
        p = subprocess.Popen(solver_cmd(run), cwd=case, stdout=fo, stderr=subprocess.STDOUT)
        pos, asked = 0, None
        while p.poll() is None:
            time.sleep(POLL_S)
            with open(log, errors="replace") as fh:
                fh.seek(pos)
                chunk = fh.read()
                pos = fh.tell()
            st.feed(chunk)
            if asked is None:
                if tail is not None:
                    tail.poll()
                if _STOP["flag"]:
                    reason = "interrupted"
                elif (tail is not None and st.final_res and max(st.final_res.values()) <= run.res_ctrl
                      and tail.spread() <= study.SETTLE_TOL):
                    reason = "settled"                   # A5: residual target met, L/D settled
                elif st.stalled():
                    reason = "stalled"
                elif st.nonfinite:
                    reason = "diverged"
                if reason:
                    foamcase.set_control(case, run, stop="writeNow", cofr=cofr, diag=diag)
                    asked = time.time()
            elif time.time() - asked > 600:
                p.terminate()
                try:
                    p.wait(60)
                except subprocess.TimeoutExpired:
                    p.kill()
        rc = p.wait()
    with open(log, errors="replace") as fh:
        fh.seek(pos)
        st.feed(fh.read())
    st.finish()
    state["wall_s"] = state.get("wall_s", 0.0) + time.time() - t0
    return rc, st, reason


def execute(run, sweep=None, index=None):
    """Run (or resume) one simulation and store its result.  Exit code: 0 result
    stored or already there, 2 error recorded, 3 interrupted (resume later), 4 gave up,
    5 running in another job (nothing done)."""
    tag = run.tag
    if not acquire_lock(tag):
        info = read_lock(tag) or {}
        print(f"RUNNING ELSEWHERE {tag}: job {_slurm_name(info.get('slurm') or {})} on "
              f"{info.get('host')} since {info.get('time')}; nothing done", flush=True)
        return 5
    try:
        return _execute(run, sweep, index)
    finally:
        release_lock(tag)


def _execute(run, sweep=None, index=None):
    tag = run.tag
    set_aside_stale(tag)
    if os.path.exists(result_path(tag)):
        print(f"done already: {tag}", flush=True)
        return 0
    if len(error_paths(tag)) >= MAX_ERRORS:
        print(f"GAVE UP after {MAX_ERRORS} errors: {tag} (see {RUNS}/{tag}.error-*.json)", flush=True)
        return 4
    if shutil.which("simpleFoam") is None:
        write_error(run, "environment", "simpleFoam not on PATH: OpenFOAM module not loaded")
        return 2
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", run.np))
    if cpus < run.np:
        write_error(run, "environment", f"run needs {run.np} cores, job has {cpus}")
        return 2
    os.makedirs(RUNS, exist_ok=True)
    case = os.path.join(CASES, tag)
    spath = os.path.join(case, "STATE.json")
    signal.signal(signal.SIGUSR1, _on_usr1)
    signal.signal(signal.SIGTERM, _on_usr1)

    # ---- set up (mesh), unless resuming
    if os.path.exists(spath):
        state = json.load(open(spath))
        print(f"resuming {tag} from iteration {state.get('resume_from', 0)} "
              f"(segment {state['segments']})", flush=True)
    else:
        try:
            grid, chk = write_case(run, case, diag=sweep in DIAG_SWEEPS)
            cm = check_mesh(case, grid)
            if cm["failures"]:
                raise RuntimeError("checkMesh: " + "; ".join(cm["failures"]))
            seq = sequence_start(run, case) if study.sequence_run(run) is not None else None
            if run.np > 1:
                sh(["decomposePar", "-force"], case, "log.decomposePar")
        except SequenceInterrupted:
            shutil.rmtree(case, ignore_errors=True)                # starts afresh next time
            print(f"INTERRUPTED {tag} during its first run; starts afresh when resubmitted", flush=True)
            return 3
        except Exception as exc:  # mesh or set-up failure: an error, never a result
            write_error(run, "mesh", str(exc), case)
            return 2
        qc = grid.foil.quarter_chord()
        state = dict(tag=tag, run=run.to_dict(), segments=0, resume_from=0,
                     cpu_s=seq["cpu_s"] if seq else 0.0, wall_s=seq["wall_s"] if seq else 0.0,
                     cells=cm["cells"], generator=grid.info, check=chk, checkMesh=cm,
                     cofr=[float(qc[0]), float(qc[1])], started=now(), sweeps=[],
                     diag=sweep in DIAG_SWEEPS, sequence=seq, long=sweep in LONG_SWEEPS)
        _atomic_json(spath, state)
    if sweep and [sweep, index] not in state["sweeps"]:
        state["sweeps"].append([sweep, index])

    # ---- solve, segment by segment
    rc, st, reason = run_segment(run, case, state, state["cofr"])
    seg_log = os.path.join(case, f"log.simpleFoam.{state['segments']}")
    if rc != 0:
        kind = foamio.log_failure_kind(seg_log)
        if kind == "numerical" and st.started:
            reason = "diverged"                             # a property of the setting: a result
        elif reason in ("stalled", "interrupted", "settled") and kind == "io":
            rc = 0                                          # it stopped where we asked it to
        else:
            # a set-up, launch or file failure, not a property of the setting: an error
            where = f"after iteration {st.it}" if st.started else "before its first iteration"
            write_error(run, "solver", f"simpleFoam exited with {rc} {where} ({kind or 'no message'})", case)
            if state["segments"] == 0 and not st.started:
                shutil.rmtree(case, ignore_errors=True)       # the next attempt starts afresh
            elif st.started:
                latest = foamio.latest_time(case if run.np == 1 else os.path.join(case, "processor0"))
                state["resume_from"] = int(float(latest)) if latest else 0
                state["cpu_s"] += foamio.read_log(seg_log, upto=state["resume_from"]).exec_s
                state["segments"] += 1
                _atomic_json(spath, state)                  # the next attempt resumes
            return 2
    if reason == "interrupted":
        latest = foamio.latest_time(case if run.np == 1 else os.path.join(case, "processor0"))
        resume = int(float(latest)) if latest else 0
        state["cpu_s"] += foamio.read_log(seg_log, upto=resume).exec_s
        state["segments"] += 1
        state["resume_from"] = resume
        _atomic_json(spath, state)
        print(f"INTERRUPTED {tag} at iteration {st.it}; resumes from {resume}", flush=True)
        return 3
    state["cpu_s"] += st.exec_s
    state["segments"] += 1
    _atomic_json(spath, state)
    return finish(run, case, state, st, reason, rc)


def finish(run, case, state, st, reason, rc):
    """Read the forces, y+ and costs, write the result JSON and the compact log, clean up."""
    tag = run.tag
    errors = []
    iters = st.it
    converged = bool(st.converged) or reason == "settled"
    if reason == "diverged" or rc != 0:
        status, stop = "failed", "diverged"
    elif reason == "settled":
        status, stop = "ok", "settled"                 # A5: residual target met and L/D settled
    elif converged:
        status, stop = "ok", "converged"
    elif reason == "stalled":
        status, stop = "ok", "stalled"
    else:
        status, stop = "ok", "max_iter" if iters >= run.max_iter else "unknown"
    res = dict(schema=1, mesh_version=MESH_VERSION, setup_version=study.SETUP_VERSION, tag=tag, status=status,
               converged=converged and status == "ok",
               stop_reason=stop, run=run.to_dict(), numerics=run.num,
               design=dict(m=run.m, t=run.t, p=study.P_CAMBER, alpha_deg=study.ALPHA_DEG,
                           trailing_edge="closed (a4 = -0.1036)"),
               setting=run.setting(), np=run.np, iterations=iters, max_iter=run.max_iter,
               residuals_final=st.final_res, cells=state["cells"],
               mesh=dict(generator=state["generator"], check=state["check"], checkMesh=state["checkMesh"]),
               sweeps=state["sweeps"], segments=state["segments"])
    # forces
    try:
        t, c = foamio.read_coefficients(case)
        cd, cl = c["Cd"], c["Cl"]
        cm = c.get("CmPitch", np.full_like(cd, np.nan))
        ld = cl / cd
        n_win = max(50, len(t) // 10)
        if state.get("long"):
            hist = os.path.join(RUNS, "history", f"{tag}.csv")
            os.makedirs(os.path.dirname(hist), exist_ok=True)
            np.savetxt(hist + ".tmp", np.column_stack([t, cl, cd, cm]), delimiter=",",
                       header="iteration,cl,cd,cm", comments="", fmt=["%d", "%.10g", "%.10g", "%.10g"])
            os.replace(hist + ".tmp", hist)
            res["history_file"] = os.path.relpath(hist, RUNS)
        res.update(cl=float(cl[-1]), cd=float(cd[-1]), cm=float(cm[-1]), ld=float(ld[-1]),
                   ld_window=dict(n=int(min(n_win, len(t))), min=float(ld[-n_win:].min()),
                                  max=float(ld[-n_win:].max())),
                   history=dict(it=foamio.decimate(t.astype(int).tolist()),
                                cl=foamio.decimate(cl.tolist()), cd=foamio.decimate(cd.tolist())))
    except Exception as exc:
        errors.append(f"forces: {exc}")
        res.update(cl=math.nan, cd=math.nan, cm=math.nan, ld=math.nan)
        if status == "ok":
            res["status"], res["stop_reason"], res["converged"] = "failed", "no_forces", False
    res["residual_history"] = foamio.decimate(st.history)
    res["residual_fields"] = foamio.residual_fields(st.field_history)
    # y+ and the first-cell height on the wall, from the final fields
    if res["status"] == "ok":
        try:
            if run.np > 1:
                sh(["reconstructPar", "-latestTime"], case, "log.reconstructPar")
            # runs the function objects of controlDict at the last time: yPlus, wallShearStress, C
            sh(["simpleFoam", "-postProcess", "-latestTime"], case, "log.postProcess")
            tdir = os.path.join(case, foamio.latest_time(case))
            yp = foamio.patch_scalars(os.path.join(tdir, "yPlus"), "airfoil")
            tw = foamio.patch_vectors(os.path.join(tdir, "wallShearStress"), "airfoil")
            cf = foamio.patch_vectors(os.path.join(tdir, "C"), "airfoil")
            utau = np.sqrt(np.linalg.norm(tw, axis=1))
            yc = yp * study.NU / np.maximum(utau, 1e-300)
            order = np.argsort(np.arctan2(cf[:, 1] - 0.0, cf[:, 0] - 0.5))
            res["yplus"] = dict(n_faces=int(len(yp)), min=float(yp.min()), mean=float(yp.mean()),
                                p95=float(np.percentile(yp, 95)), max=float(yp.max()),
                                frac_gt_1=float(np.mean(yp > 1.0)))
            res["first_cell"] = dict(min=float(2 * yc.min()), median=float(2 * np.median(yc)),
                                     max=float(2 * yc.max()), target=run.dy1)
            sel = foamio.decimate(order.tolist(), 120)
            res["wall_profile"] = dict(x=cf[sel, 0].tolist(), y=cf[sel, 1].tolist(),
                                       yplus=yp[sel].tolist(), cf=(2 * np.linalg.norm(tw[sel], axis=1)).tolist())
        except Exception as exc:
            errors.append(f"yplus: {exc}")
    if state.get("diag") and run.np == 1:
        try:
            res["fields"] = field_diagnostics(run, case)
        except Exception as exc:
            errors.append(f"fields: {exc}")
    # costs and provenance
    log0 = os.path.join(case, "log.simpleFoam.0")
    head = open(log0, errors="replace").read(20000) if os.path.exists(log0) else ""
    seq = state.get("sequence")
    # A4: the first run of a sequenced start is part of the run: its cost is added
    res["cost"] = dict(cells_iters=float(state["cells"]) * iters + (seq["cells_iters"] if seq else 0.0),
                       cpu_s=state["cpu_s"], wall_s=state["wall_s"],
                       core_s=(state["cpu_s"] - (seq["cpu_s"] if seq else 0.0)) * run.np + (seq["cpu_s"] if seq else 0.0))
    if seq:
        res["sequence"] = seq
    res["openfoam"] = dict(version=foamio.foam_version(head) or os.environ.get("WM_PROJECT_VERSION"),
                           env_version=os.environ.get("WM_PROJECT_VERSION"),
                           module=os.environ.get("AIRFOIL_OF_MODULE"),
                           decay_control_echoed="decayControl" in head)
    res.update(host=socket.gethostname(), slurm=slurm_ids(), code=code_hash(),
               started=state["started"], finished=now(), errors=errors)
    _atomic_json(result_path(tag), res)
    # compact log: set-up logs and the head and tail of every solver segment
    with open(os.path.join(RUNS, f"{tag}.log.txt"), "w") as fo:
        for name in sorted(os.listdir(case)):
            if not name.startswith("log."):
                continue
            lines = open(os.path.join(case, name), errors="replace").readlines()
            fo.write(f"==> {name} ({len(lines)} lines)\n")
            fo.writelines(lines if len(lines) <= 400 else lines[:150] + ["   ...\n"] + lines[-250:])
    shutil.rmtree(case, ignore_errors=True)
    print(f"{res['status'].upper()} {tag}: L/D {res['ld']:.6g} cl {res['cl']:.6g} cd {res['cd']:.6g} "
          f"iterations {iters} ({res['stop_reason']}), {state['cells']} cells, "
          f"{state['cpu_s'] / 60:.1f} cpu-min", flush=True)
    return 0


def mesh_only(run):
    """The check job's mesh check: generator check, polyMesh, checkMesh; exit 0 if it passes."""
    case = os.path.join(CASES, "meshcheck_" + run.tag)
    out = dict(tag=run.tag, mesh_version=MESH_VERSION, setting=run.setting(), design=[run.m, run.t],
               host=socket.gethostname(), time=now())
    try:
        t0 = time.time()
        grid, chk = write_case(run, case)
        out.update(cells=grid.n_cells(), generator=grid.info, check=chk, write_s=time.time() - t0)
        cm = check_mesh(case, grid)
        out["checkMesh"] = cm
        out["passed"] = not cm["failures"]
    except Exception as exc:
        out.update(passed=False, error=str(exc))
    d = os.path.join(WORK, "meshcheck")
    _atomic_json(os.path.join(d, f"{run.tag}.json"), out)
    log = os.path.join(case, "log.checkMesh")
    if os.path.exists(log):
        shutil.copy(log, os.path.join(d, f"{run.tag}.checkMesh.txt"))
    shutil.rmtree(case, ignore_errors=True)
    cm = out.get("checkMesh", {})
    print(f"MESH {'PASSED' if out['passed'] else 'FAILED'} {run.tag}: {out.get('cells')} cells, checkMesh "
          f"non-orthogonality {cm.get('max_non_orth')}, skewness {cm.get('max_skewness')}, "
          f"failures {cm.get('failures', out.get('error'))}", flush=True)
    return 0 if out["passed"] else 2


# ---------------------------------------------------------------- listing
def resource_class(run):
    """(time limit, memory) for the Slurm array of a run of this size."""
    if run.res_ctrl <= study.SETTLE_RES_CTRL:
        # the settling runs: 40,000 iterations of the reference mesh, about 20 hours at the
        # reference's measured 2.1 us per cell-iteration and 43 at 4.6; past the limit the
        # USR1 checkpoint stops them cleanly and a resubmission resumes them
        t, mem = "2-12:00:00", 8
    elif run.n_surf <= 128:
        t, mem = "03:00:00", 2
    elif run.n_surf <= 181:
        t, mem = "06:00:00", 3
    elif run.n_surf <= 256:
        t, mem = "12:00:00", 4
    elif run.n_surf <= 362:
        t, mem = "1-12:00:00", 8
    else:
        t, mem = "3-00:00:00", 12
    return t, f"{mem + run.np - 1}G"


def compress(idx):
    """[0, 1, 2, 5, 7, 8] -> '0-2,5,7-8' (a Slurm --array list)."""
    out, idx = [], sorted(idx)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
            j += 1
        out.append(f"{idx[i]}" if i == j else f"{idx[i]}-{idx[j]}")
        i = j + 1
    return ",".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep")
    ap.add_argument("--index", type=int)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--json", action="store_true", help="with --list: one JSON object per line")
    ap.add_argument("--count", action="store_true")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--groups", action="store_true", help="unfinished indices by resource class")
    ap.add_argument("--dry-run", action="store_true", help="write the case directory and stop")
    ap.add_argument("--mesh-only", action="store_true",
                    help="write the mesh, run checkMesh, record $AIRFOIL_WORK/meshcheck/<tag>.json, clean up")
    ap.add_argument("--case-dir", help="with --dry-run: where to write the case")
    ap.add_argument("--run-setting", nargs="+", metavar="KEY=VALUE",
                    help="one run outside any sweep: m, t, n_surf, dy1, radius, res_ctrl, gamg_tol [, np]")
    ap.add_argument("--design-from-config", action="store_true",
                    help="with --run-setting: take m and t from config/stage1.json (the hardest design)")
    ap.add_argument("--print", action="store_true", help="with --run-setting: print L/D and cost for errbudget")
    a = ap.parse_args(argv)

    if a.run_setting:
        kv = dict(x.split("=", 1) for x in a.run_setting)
        if a.design_from_config:
            kv["m"], kv["t"] = study.hardest_design()
        run = study.run_at((float(kv["m"]), float(kv["t"])), kv, np=int(kv.get("np", 1)))
        runs, idx = [run], 0
    else:
        if not a.sweep:
            ap.error("give --sweep NAME or --run-setting")
        try:
            runs = study.sweep(a.sweep)
        except (RuntimeError, KeyError) as exc:
            print(f"{a.sweep}: {exc}", file=sys.stderr)
            return 1
        idx = a.index

    if a.count:
        print(len(runs))
        return 0
    if a.list:
        from af import costmodel
        for i, r in enumerate(runs):
            c = costmodel.cells(r.n_surf, r.dy1, r.radius)
            h = costmodel.predicted_core_seconds(r) / 3600
            row = dict(index=i, tag=r.tag, np=r.np, cells=c, predicted_core_h=round(h, 3),
                       status=status_of(r), time=resource_class(r)[0], mem=resource_class(r)[1])
            print(json.dumps(row) if a.json else
                  f"{i:4d}  {r.tag:<58s} np={r.np}  cells={c:>8d}  ~{h:6.2f} core-h  {row['status']}")
        return 0
    if a.status:
        counts = {}
        for r in runs:
            s = status_of(r)
            counts[s] = counts.get(s, 0) + 1
        print(f"{a.sweep}: {len(runs)} runs: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
        return 0
    if a.groups:
        groups = {}
        for i, r in enumerate(runs):
            if status_of(r) in ("done", "failed", "gave-up", "running"):
                continue
            key = (resource_class(r)[0], resource_class(r)[1], r.np)
            groups.setdefault(key, []).append(i)
        for (t, mem, np_), ii in sorted(groups.items(), key=lambda kv: kv[1][0]):
            print(f"{t} {mem} {np_} {compress(ii)}")
        return 0

    if idx is None or not (0 <= idx < len(runs)):
        ap.error(f"--index must be in 0..{len(runs) - 1}")
    run = runs[idx]
    if a.mesh_only:
        return mesh_only(run)
    if a.dry_run:
        case = a.case_dir or os.path.join(WORK, "dryrun", f"{a.sweep or 'single'}_{idx}_{run.tag}")
        grid, chk = write_case(run, case, diag=a.sweep in DIAG_SWEEPS)
        print(f"dry run: {run.tag}\n  case written to {case}\n  {grid.n_cells()} cells; "
              f"generator check: {json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in chk.items()})}")
        return 0
    code = execute(run, a.sweep, idx)
    if a.print:
        p = result_path(run.tag)
        if not os.path.exists(p):
            print("RUN FAILED (no result)")
            return 1
        r = json.load(open(p))
        if r["status"] != "ok" or not math.isfinite(r.get("ld", math.nan)):
            print("RUN FAILED")
            return 1
        print(f"L/D = {r['ld']!r}\ncl = {r['cl']!r}\ncd = {r['cd']!r}\ncost = {r['cost']['cells_iters']!r}")
        print("converged = yes" if r["converged"] else "NOT CONVERGED")
        return 0
    return code


if __name__ == "__main__":
    sys.exit(main())
