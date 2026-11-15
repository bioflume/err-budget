"""Readers for what OpenFOAM writes: solver logs, forceCoeffs output, checkMesh, patch fields.

The forceCoeffs reader follows hamilton/common/ofmodern.py of the article's
repository (columns found by name, every restart segment read, the last value of a
repeated time kept); v2206 and later write coefficient.dat with '# Time Cd ... Cl'.
"""
from __future__ import annotations

import glob
import math
import os
import re

import numpy as np

NUM_DIR = re.compile(r"[0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)?")

_TIME = re.compile(r"^Time = ([0-9.eE+-]+)s?\s*$")
_SOLVE = re.compile(r"Solving for (\w+), Initial residual = ([0-9.eE+-]+|nan|-?inf)", re.I)
_EXEC = re.compile(r"^ExecutionTime = ([0-9.eE+-]+) s\s+ClockTime = ([0-9.eE+-]+) s")


class LogState:
    """Incremental reader of a simpleFoam log: iteration, residuals, stall detection."""

    def __init__(self, stall_iter: int, first_iter: int = 0):
        self.it = first_iter
        self.res = {}              # residuals of the current iteration (first solve of each field)
        self.history = []          # (iteration, max initial residual)
        self.field_history = []    # (iteration, {field: initial residual})
        self.best = math.inf
        self.best_it = first_iter
        self.stall_iter = stall_iter
        self.exec_s = 0.0
        self.clock_s = 0.0
        self.converged = False
        self.nonfinite = False
        self.started = False
        self.final_res = {}
        self._buf = ""

    def feed(self, text: str):
        self._buf += text
        lines = self._buf.split("\n")
        self._buf = lines.pop()
        for line in lines:
            self.line(line)

    def line(self, line: str):
        m = _TIME.match(line.strip())
        if m:
            self._close()
            self.it = int(round(float(m.group(1))))
            self.started = True
            self.res = {}
            return
        m = _SOLVE.search(line)
        if m and m.group(1) not in self.res:
            try:
                v = float(m.group(2))
            except ValueError:
                v = math.nan
            self.res[m.group(1)] = v
            return
        m = _EXEC.match(line.strip())
        if m:
            self.exec_s, self.clock_s = float(m.group(1)), float(m.group(2))
            return
        if "solution converged in" in line:
            self.converged = True

    def _close(self):
        if not self.res:
            return
        vals = list(self.res.values())
        if not all(math.isfinite(v) for v in vals):
            self.nonfinite = True
        r = max(vals)
        self.final_res = dict(self.res)
        self.history.append((self.it, r))
        self.field_history.append((self.it, dict(self.res)))
        if r < 0.95 * self.best:
            self.best, self.best_it = r, self.it

    def finish(self):
        self._close()
        self.res = {}

    def stalled(self) -> bool:
        return self.started and self.it - self.best_it >= self.stall_iter


def read_log(path, stall_iter=10**9, upto=None):
    """A whole log: LogState after the last iteration at or before `upto` (all if None)."""
    st = LogState(stall_iter)
    if not os.path.exists(path):
        return st
    with open(path, errors="replace") as fh:
        for line in fh:
            m = _TIME.match(line.strip())
            if m and upto is not None and float(m.group(1)) > upto + 1e-9:
                break
            st.line(line)
    st.finish()
    return st


def log_failure_kind(path) -> str | None:
    """How a solver run that exited with an error failed:
      'numerical'  a floating-point exception, a segmentation fault, or a fatal error of
                   the solution (not of reading a file) after the first iteration;
      'io'         a fatal error reading or writing a file (e.g. a dictionary caught
                   half-written), which says nothing about the setting;
      'setup'      any fatal error before the first iteration;
      None         no failure visible in the log."""
    if not os.path.exists(path):
        return "setup"
    txt = open(path, errors="replace").read()
    started = re.search(r"^Time = ", txt, flags=re.M) is not None
    crash = ("sigFpe" in txt or "Floating point exception" in txt or "Segmentation fault" in txt
             or "sigSegv" in txt)
    io = "FOAM FATAL IO ERROR" in txt
    fatal = "FOAM FATAL ERROR" in txt
    if not (crash or io or fatal):
        return None
    if not started:
        return "setup"
    if crash or (fatal and not io):
        return "numerical"
    return "io"


def foam_version(log_text: str) -> str | None:
    m = re.search(r"Version:\s+(\S+)", log_text)
    return m.group(1) if m else None


def read_coefficients(case, name="forceCoeffs"):
    """(iteration, {column: array}) from forceCoeffs output of every restart segment."""
    files = sorted(glob.glob(os.path.join(case, "postProcessing", name, "*", "coefficient*.dat"))
                   + glob.glob(os.path.join(case, "postProcessing", name, "*", "forceCoeffs*.dat")),
                   key=lambda p: float(os.path.basename(os.path.dirname(p))))
    cols, rows = None, []
    for f in files:
        hdr = None
        for line in open(f):
            if line.startswith("#"):
                toks = line[1:].split()
                if toks and toks[0] == "Time":
                    hdr = toks
                continue
            if hdr and line.strip():
                parts = line.split()
                if len(parts) >= len(hdr):
                    rows.append({h: float(v) for h, v in zip(hdr, parts)})
        cols = cols or hdr
    if not rows:
        raise RuntimeError(f"no forceCoeffs output in {case}/postProcessing/{name}")
    rows.sort(key=lambda r: r["Time"])
    keep = {}
    for r in rows:
        keep[r["Time"]] = r                     # restarts repeat times: keep the last
    t = np.array(sorted(keep))
    out = {k: np.array([keep[x].get(k, np.nan) for x in t]) for k in cols}
    return t, out


# a number as OpenFOAM prints it; never the full stop that checkMesh writes after some
# ("Min volume = 4.17e-08. Max volume = ...")
NUM = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"


def parse_check_mesh(text: str) -> dict:
    out = {}
    for key, pat in (("cells", r"^\s*cells:\s+(\d+)"), ("faces", r"^\s*faces:\s+(\d+)"),
                     ("points", r"^\s*points:\s+(\d+)"),
                     ("max_non_orth", r"Mesh non-orthogonality Max:\s*" + NUM),
                     ("avg_non_orth", r"Mesh non-orthogonality Max:\s*[-+0-9.eE]+?\s+average:\s*" + NUM),
                     ("max_skewness", r"Max skewness\s*=\s*" + NUM),
                     ("max_aspect", r"Max aspect ratio\s*[=:]\s*" + NUM),
                     ("min_volume", r"Min volume\s*=\s*" + NUM),
                     ("max_volume", r"Max volume\s*=\s*" + NUM)):
        m = re.search(pat, text, flags=re.M)
        if m:
            out[key] = float(m.group(1)) if key not in ("cells", "faces", "points") else int(m.group(1))
    out["mesh_ok"] = "Mesh OK." in text
    m = re.search(r"Failed (\d+) mesh checks", text)
    out["failed_checks"] = int(m.group(1)) if m else 0
    out["fatal"] = bool(re.search(r"\*\*\*(Zero or negative cell volume|Error in face pyramids|"
                                  r"Incorrectly oriented faces)", text))
    return out


# checkMesh's high-aspect-ratio check fails at an aspect ratio of 1000, a fixed threshold
# that every wall-resolved mesh (y+ ~ 1) passes; it is the one failed check allowed.
ALLOWED_CHECKMESH = ("High aspect ratio cells",)


def check_mesh_failures(text: str) -> list:
    """The failed checks of a checkMesh log that fail a run (PROTOCOL.md 14, A2): every
    '***' line except the allowed ones, and severely non-orthogonal faces."""
    out = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("***") and not any(a in s for a in ALLOWED_CHECKMESH):
            out.append(s)
        elif "severely non-orthogonal" in s:
            out.append(s)
    m = re.search(r"Failed (\d+) mesh checks", text)
    allowed = sum(1 for line in text.splitlines()
                  if line.strip().startswith("***") and any(a in line for a in ALLOWED_CHECKMESH))
    if m and int(m.group(1)) > allowed and not out:
        out.append(f"checkMesh: Failed {m.group(1)} mesh checks")
    return out


def _patch_block(text, patch):
    i = text.index("boundaryField")
    m = re.search(r"\b%s\b\s*\{" % re.escape(patch), text[i:])
    if not m:
        raise KeyError(patch)
    return text[i + m.end():]


def patch_scalars(path, patch):
    blk = _patch_block(open(path).read(), patch)
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", blk)
    if m:
        n = int(m.group(1))
        return np.array(blk[m.end():].split(")", 1)[0].split()[:n], float)
    m = re.search(r"value\s+uniform\s+([-0-9.eE+]+)", blk)
    return np.array([float(m.group(1))])


def patch_vectors(path, patch):
    blk = _patch_block(open(path).read(), patch)
    m = re.search(r"value\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", blk)
    if not m:
        raise ValueError(f"no nonuniform vector list for {patch} in {path}")
    n = int(m.group(1))
    v = re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)", blk[m.end():])[:n]
    return np.array(v, float)


def latest_time(case):
    ts = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", d) and os.path.isdir(os.path.join(case, d))]
    return max(ts, key=float) if ts else None


def decimate(x, n=200):
    """At most n roughly evenly spaced entries, always with the last one."""
    x = list(x)
    if len(x) <= n:
        return x
    idx = sorted(set(np.linspace(0, len(x) - 1, n).round().astype(int).tolist()))
    return [x[i] for i in idx]


def residual_fields(field_history, n=300):
    """{"it": [...], field: [...]} at most n entries, from LogState.field_history."""
    rows = decimate(field_history, n)
    names = sorted({k for _, r in rows for k in r})
    out = {"it": [it for it, _ in rows]}
    for k in names:
        out[k] = [r.get(k, math.nan) for _, r in rows]
    return out


def internal_field(path, n=None):
    """The internalField of an ascii volScalarField or volVectorField: (n,) or (n, 3)."""
    with open(path, errors="replace") as fh:
        text = fh.read()
    i = text.index("internalField")
    head = text[i:i + 200]
    m = re.match(r"internalField\s+uniform\s+(\([^)]*\)|[^;\s]+)\s*;", head)
    if m:
        v = np.array([float(x) for x in m.group(1).strip("()").split()])
        if n is None:
            raise ValueError(f"{path}: uniform field and no cell count")
        return np.tile(v, (n, 1)) if v.size > 1 else np.full(n, v[0])
    m = re.match(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", head)
    if not m:
        raise ValueError(f"{path}: internalField not understood")
    kind, count = m.group(1), int(m.group(2))
    start = i + m.end()
    stop = text.index("\n)", start)
    body = text[start:stop].replace("(", " ").replace(")", " ")
    vals = np.array(body.split(), dtype=float)
    if kind == "vector":
        vals = vals.reshape(-1, 3)
    if len(vals) != count:
        raise ValueError(f"{path}: {len(vals)} values, header says {count}")
    return vals


class CoeffTail:
    """A5: reads the forceCoeffs output of a running case as it grows and keeps L/D over the
    last `window` + 1 iterations, to tell when it has settled."""

    def __init__(self, case, window, name="forceCoeffs"):
        self.case, self.window, self.name = case, int(window), name
        self.pos, self.hdr = {}, {}
        self.it, self.ld = [], []

    def poll(self):
        for f in glob.glob(os.path.join(self.case, "postProcessing", self.name, "*", "coefficient*.dat")):
            with open(f, errors="replace") as fh:
                fh.seek(self.pos.get(f, 0))
                chunk = fh.read()
                end = chunk.rfind("\n")
                if end < 0:
                    continue
                self.pos[f] = self.pos.get(f, 0) + end + 1
            for line in chunk[:end].split("\n"):
                if line.startswith("#"):
                    toks = line[1:].split()
                    if toks and toks[0] == "Time":
                        self.hdr[f] = toks
                    continue
                h = self.hdr.get(f)
                parts = line.split()
                if not h or len(parts) < len(h):
                    continue
                row = dict(zip(h, parts))
                try:
                    i, cl, cd = int(round(float(row["Time"]))), float(row["Cl"]), float(row["Cd"])
                except (KeyError, ValueError):
                    continue
                if self.it and i <= self.it[-1]:        # a restart repeats iterations
                    while self.it and self.it[-1] >= i:
                        self.it.pop()
                        self.ld.pop()
                self.it.append(i)
                self.ld.append(cl / cd if cd else math.nan)
        keep = self.window + 1
        if len(self.ld) > keep:
            del self.it[:-keep], self.ld[:-keep]

    def spread(self):
        """(max - min) / |last| of L/D over the last window + 1 iterations, or inf."""
        if len(self.ld) < self.window + 1 or not all(math.isfinite(v) for v in self.ld):
            return math.inf
        return (max(self.ld) - min(self.ld)) / abs(self.ld[-1])
