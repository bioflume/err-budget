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
import subprocess

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


def run(app, case, nprocs=1, log=None, extra=()):
    """Run an OpenFOAM application, in parallel if nprocs > 1.  Returns stdout."""
    if nprocs > 1:
        _w(f"{case}/system/decomposeParDict", "dictionary", "decomposeParDict",
           decompose_dict(nprocs))
        if not glob.glob(os.path.join(case, "processor0")):
            _sh(["decomposePar", "-force"], case)
        cmd = ["mpirun", "-np", str(nprocs), app, "-parallel", *extra]
    else:
        cmd = [app, *extra]
    with open(log or os.devnull, "w") as fo:
        r = subprocess.run(of_cmd(cmd), cwd=case, env=ENV, stdout=fo,
                           stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        tail = open(log).read()[-1500:] if log else ""
        raise RuntimeError(f"{' '.join(cmd)} failed in {case}:\n{tail}\n{r.stderr[-1500:]}")


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
