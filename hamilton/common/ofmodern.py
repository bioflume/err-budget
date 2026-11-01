"""Helpers for current OpenFOAM releases (function objects, parallel runs).

The Ubuntu package of OpenFOAM v1912 used for the article cannot construct
function objects, so the article's drivers compute forces from written fields.
On Hamilton a current release is used and the forces come from OpenFOAM's own
forceCoeffs function object; the field-based post-processing is kept as a
check (J1).
"""
from __future__ import annotations

import glob
import os
import re
import shlex
import signal
import subprocess
import tempfile
import time

import numpy as np

from openfoam_cylinder import ENV, of_cmd, _w


def force_coeffs_fo(name="fc", patch="cylinder", magU=1.0, lref=0.1, aref=0.001,
                    cofr=(0.2, 0.2, 0.005), every=1):
    """A forceCoeffs function object writing every `every` time steps."""
    return f"""    {name}
    {{
        type            forceCoeffs;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   {every};
        log             false;
        patches         ({patch});
        rho             rhoInf;
        rhoInf          1;
        CofR            ({cofr[0]} {cofr[1]} {cofr[2]});
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
        magUInf         {magU};
        lRef            {lref};
        Aref            {aref};
    }}
"""


def control_dict(start, end, dt, write_interval, functions="", app="pimpleFoam",
                 write_control="runTime"):
    fn = f"functions\n{{\n{functions}}}\n" if functions else ""
    return f"""application {app};
startFrom latestTime; startTime {start:.10g}; stopAt endTime; endTime {end:.10g};
deltaT {dt:.10g}; writeControl {write_control}; writeInterval {write_interval:.10g};
purgeWrite 0; writeFormat ascii; writePrecision 16; writeCompression off;
timeFormat general; timePrecision 10; runTimeModifiable false; adjustTimeStep no;
{fn}"""


def decompose_dict(n):
    """scotch (default); ERRBUDGET_OF_DECOMP=simple for builds without scotch."""
    m = os.environ.get("ERRBUDGET_OF_DECOMP", "scotch")
    if m == "simple":
        return (f"numberOfSubdomains {n};\nmethod simple;\n"
                f"simpleCoeffs {{ n ({n} 1 1); delta 0.001; }}\n")
    return f"numberOfSubdomains {n};\nmethod {m};\n"


class TooSlow(RuntimeError):
    """A run stopped by the watchdog of run(): the machine, not the run, is at fault."""


_ITERS = re.compile(r"No Iterations (\d+)")


def _watch(p, log, limit, window, poll=60.0):
    """Wait for the process p.  Every `poll` seconds, count the linear-solver work that the
    log reports, each solve as one unit plus its iterations (the "No Iterations N" of every
    "Solving for" line), and compare it with the time of the last `window` seconds.
    Return None when p ends, or a description if a unit took more than `limit` seconds.

    On Hamilton in October 2026, the J2 spin-up of n_arc 192 (32 ranks) ran at about
    4 s per linear-solver iteration, against 2-5 ms on the coarser meshes, and the J2
    reference at n_arc 128 and dt 2.5e-4 (16 ranks) took about 1,250 s per time step,
    against 0.9 s at dt 1.25e-4 on another node; both used up their time limits (seven
    days for the first).  Counting solver work separates such a slowdown from a hard
    numerical problem, whose solves take many iterations at the normal cost each; a
    solve that needs no iteration still counts, so that loose tolerances do not look
    slow."""
    t0, pos, rest, total = time.time(), 0, "", 0
    hist = [(t0, 0)]
    while True:
        try:
            p.wait(timeout=poll)
            return None
        except subprocess.TimeoutExpired:
            pass
        try:
            with open(log, "rb") as f:
                f.seek(pos)
                chunk = f.read().decode(errors="replace")
                pos = f.tell()
        except OSError:
            chunk = ""
        text = rest + chunk
        cut = text.rfind("\n") + 1
        rest = text[cut:]
        total += sum(1 + int(m) for m in _ITERS.findall(text[:cut]))
        now = time.time()
        hist.append((now, total))
        old = [h for h in hist if now - h[0] >= window]
        if not old:
            continue
        t_old, n_old = old[-1]
        hist = [h for h in hist if h[0] >= t_old]
        n = total - n_old
        if now - t_old > limit * max(n, 1):
            return (f"{n} linear solves and iterations in the last {(now - t_old) / 60:.1f} "
                    f"minutes, {(now - t_old) / max(n, 1):.3g} s each (limit {limit:g} s)")


def _kill(p):
    """Stop p and everything it started (it leads its own process group)."""
    for sig, wait in ((signal.SIGTERM, 60), (signal.SIGKILL, 30)):
        try:
            os.killpg(p.pid, sig)
        except (ProcessLookupError, PermissionError):
            return
        try:
            p.wait(timeout=wait)
            return
        except subprocess.TimeoutExpired:
            pass


def run(app, case, nprocs=1, log=None, extra=()):
    """Run an OpenFOAM application, in parallel if nprocs > 1.

    ERRBUDGET_OF_MPIRUN_OPTS: extra mpirun options (hamilton/config.sh sets
    --bind-to none).  ERRBUDGET_OF_MAX_S_PER_ITER > 0 (with a log): stop the run and
    raise TooSlow when, over the last ERRBUDGET_OF_WATCH_WINDOW seconds (default 900),
    its linear solves and their iterations took longer than that each (_watch)."""
    if nprocs > 1:
        _w(f"{case}/system/decomposeParDict", "dictionary", "decomposeParDict",
           decompose_dict(nprocs))
        if not glob.glob(os.path.join(case, "processor0")):
            _sh(["decomposePar", "-force"], case)
        opts = shlex.split(os.environ.get("ERRBUDGET_OF_MPIRUN_OPTS", ""))
        cmd = ["mpirun", *opts, "-np", str(nprocs), app, "-parallel", *extra]
    else:
        cmd = [app, *extra]
    limit = float(os.environ.get("ERRBUDGET_OF_MAX_S_PER_ITER", "0") or 0)
    window = float(os.environ.get("ERRBUDGET_OF_WATCH_WINDOW", "900") or 900)
    slow = None
    with open(log or os.devnull, "w") as fo, tempfile.TemporaryFile("w+") as fe:
        p = subprocess.Popen(of_cmd(cmd), cwd=case, env=ENV, stdout=fo, stderr=fe, text=True,
                             start_new_session=True)
        try:
            if limit > 0 and log:
                slow = _watch(p, log, limit, window, poll=min(60.0, max(1.0, window / 5)))
                if slow:
                    _kill(p)
            else:
                p.wait()
        except BaseException:
            _kill(p)
            raise
        fe.seek(0)
        err = fe.read()[-1500:]
    if slow:
        raise TooSlow(f"{' '.join(cmd)} in {case} on {os.uname().nodename} was stopped: {slow}. "
                      "This is a slowdown of the machine (CPU binding or a busy node), not a "
                      "property of the run; the task is retried when the job is resubmitted.")
    if p.returncode != 0:
        tail = open(log, errors="replace").read()[-1500:] if log else ""
        raise RuntimeError(f"{' '.join(cmd)} failed in {case}:\n{tail}\n{err}")


def _sh(cmd, cwd):
    r = subprocess.run(of_cmd(cmd), cwd=cwd, env=ENV, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed:\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    return r.stdout


def latest_time(case, parallel=False):
    d = os.path.join(case, "processor0") if parallel else case
    ts = [t for t in os.listdir(d) if re.fullmatch(r"[\d.eE+-]+", t)]
    return max(ts, key=float) if ts else None


def read_coefficients(case, name="fc"):
    """(t, Cd, Cl) from forceCoeffs output, all restart segments, by column name.

    Current releases write postProcessing/<name>/<t0>/coefficient.dat with a header
    line '# Time Cd Cd(f) Cd(r) Cl ...'; older ones forceCoeffs.dat with
    '# Time Cm Cd Cl ...'.  Columns are found by name."""
    files = sorted(glob.glob(os.path.join(case, "postProcessing", name, "*", "coefficient*.dat"))
                   + glob.glob(os.path.join(case, "postProcessing", name, "*", "forceCoeffs*.dat")),
                   key=lambda p: float(os.path.basename(os.path.dirname(p))))
    T, CD, CL = [], [], []
    for f in files:
        hdr, rows = None, []
        for line in open(f):
            if line.startswith("#"):
                toks = line[1:].split()
                if toks and toks[0] == "Time":
                    hdr = toks
                continue
            if line.strip():
                rows.append(line.split())
        if hdr is None or not rows:
            continue
        icd, icl = hdr.index("Cd"), hdr.index("Cl")
        a = np.array([[float(r[0]), float(r[icd]), float(r[icl])] for r in rows
                      if len(r) >= len(hdr)])
        T.append(a[:, 0]); CD.append(a[:, 1]); CL.append(a[:, 2])
    if not T:
        raise RuntimeError(f"no forceCoeffs output in {case}/postProcessing/{name}")
    t, cd, cl = np.concatenate(T), np.concatenate(CD), np.concatenate(CL)
    # restarts repeat a few times: keep the last value written for each time
    order = np.argsort(t, kind="stable")
    t, cd, cl = t[order], cd[order], cl[order]
    keep = np.r_[t[1:] - t[:-1] > 1e-12, True]
    return t[keep], cd[keep], cl[keep]


def log_steps(log):
    """Per time step: time, cumulative ExecutionTime, outer iterations, and
    whether PIMPLE reported it did not converge."""
    out, cur, its, nc = [], None, 0, False
    if not os.path.exists(log):
        return out
    for line in open(log, errors="replace"):
        if line.startswith("Time = "):
            cur, its, nc = float(line.split()[2].rstrip("s")), 0, False
        elif line.startswith("PIMPLE: iteration"):
            its += 1
        elif line.startswith("PIMPLE: not converged"):
            nc = True
        elif line.startswith("ExecutionTime =") and cur is not None:
            out.append((cur, float(line.split()[2]), its, nc))
    return out
