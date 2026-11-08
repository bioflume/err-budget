"""Run a solver that is a shell command from a job file.

The engine only needs a function `solver(**params) -> (observables, cost)`.  This
module builds one around a command line: for every setting it creates a run
directory, renders `{{placeholders}}` in the command and in template files,
runs the command, reads the observable (and optionally the cost) from its
output, decides whether the run failed or did not converge, and records
everything in `runs.jsonl` so that a stopped tuning resumes without repeating
runs.

Nothing here is specific to a solver.  OpenFOAM, LAMMPS, GROMACS, a Fortran
binary or `python my_solver.py` are all "a command that reads parameters and
writes a number".
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from typing import Any, Callable

from .api import RunFailed

PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?::([^}]*))?\s*\}\}")
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
BUILTIN_PLACEHOLDERS = ("run_dir", "run_id", "job_dir")


class JobError(ValueError):
    """A job file that cannot be used as written; the message says what to fix."""


def _as_list(x) -> list:
    return [] if x is None else [x] if isinstance(x, str) else list(x)


def _to_float(s: Any) -> float:
    """float() that also accepts Fortran exponents (1.5D-03) and stray brackets."""
    if isinstance(s, (int, float)):
        return float(s)
    if s is None or isinstance(s, (list, dict)):
        raise ValueError(f"not a number: {s!r}")
    t = str(s).strip().strip("()[]{},;")
    return float(t.replace("D", "e").replace("d", "e"))


# ----------------------------------------------------------------- rendering
def format_value(v: Any, spec: str | None = None) -> str:
    """How a parameter value appears in a command or a template file.

    Integers as integers, floats in the shortest round-trip form (1e-06, 0.05),
    strings as they are.  A format spec after a colon ({{dt:.3e}}) is applied
    with str.format.
    """
    if spec:
        return format(v, spec)
    if isinstance(v, bool):
        return str(v).lower()
    if isinstance(v, float):
        return repr(v)
    return str(v)


def render(text: str, values: dict, where: str = "template") -> str:
    """Replace every {{name}} (or {{name:fmt}}) in `text`; unknown names are an error."""
    def sub(m):
        name, fmt = m.group(1), m.group(2)
        if name not in values:
            raise JobError(f"{where}: unknown placeholder {{{{{name}}}}}; known names: "
                           f"{', '.join(sorted(values))}")
        try:
            return format_value(values[name], fmt)
        except (ValueError, TypeError) as exc:
            raise JobError(f"{where}: cannot format {{{{{name}:{fmt}}}}} with value "
                           f"{values[name]!r} ({exc})") from None
    return PLACEHOLDER.sub(sub, text)


def placeholders(text: str) -> list[str]:
    return [m.group(1) for m in PLACEHOLDER.finditer(text)]


# ----------------------------------------------------------------- extraction
def _json_path(obj: Any, key: str) -> Any:
    """Walk 'a.b[2].c' through nested dicts and lists."""
    for part in re.findall(r"[^.\[\]]+|\[\d+\]", key):
        if part.startswith("["):
            obj = obj[int(part[1:-1])]
        elif isinstance(obj, dict):
            obj = obj[part]
        else:
            raise KeyError(part)
    return obj


def parse_method(spec: dict) -> str:
    p = spec.get("parse")
    if p is None:
        p = "regex" if "regex" in spec else "json" if "key" in spec else \
            "column" if "column" in spec else "last_number"
    return str(p)


def extract_number(spec: dict, text: str, what: str = "observable") -> float:
    """Read one number out of `text` as the spec says.

    spec["parse"]:
      "regex"        spec["regex"] with one capture group; the last match (or the
                     first with spec["which"] = "first")
      "json"         json.loads(text) and spec["key"] as a path ("forces.drag", "a[0].b");
                     if the text has other lines around the JSON, the last line that is
                     JSON is used
      "column"       the last non-comment line split on whitespace, commas or brackets,
                     column spec["column"] (1-based); spec["comment"] holds the
                     characters that start a comment line (default '#'; '#@' for xvg)
      "last_number"  the last number in the text (default when nothing else is given)
    Fortran exponents (1.5D-03) are understood everywhere.
    """
    parse = parse_method(spec)
    try:
        if parse == "regex":
            ms = re.findall(spec["regex"], text, re.MULTILINE)
            if not ms:
                raise JobError(f"{what}: the pattern {spec['regex']!r} matched nothing")
            m = ms[0] if spec.get("which") == "first" else ms[-1]
            if isinstance(m, tuple):
                m = m[0]
            return _to_float(m)
        if parse == "json":
            try:
                obj = json.loads(text)
            except json.JSONDecodeError:
                obj = None
                for line in reversed(text.splitlines()):
                    line = line.strip()
                    if line[:1] in ("{", "["):
                        try:
                            obj = json.loads(line)
                            break
                        except json.JSONDecodeError:
                            continue
                if obj is None:
                    raise
            return _to_float(_json_path(obj, spec["key"]))
        if parse == "column":
            comment = str(spec.get("comment", "#"))
            lines = [ln for ln in text.splitlines() if ln.strip() and ln.lstrip()[0] not in comment]
            if not lines:
                raise JobError(f"{what}: no data lines")
            row = [c for c in re.split(r"[,\s()\[\]]+", lines[-1].strip()) if c]
            return _to_float(row[int(spec["column"]) - 1])
        if parse == "last_number":
            ms = NUMBER.findall(text)
            if not ms:
                raise JobError(f"{what}: no number found in the output")
            return _to_float(ms[-1])
    except JobError:
        raise
    except (KeyError, IndexError, ValueError, TypeError, json.JSONDecodeError) as exc:
        raise JobError(f"{what}: could not read a number ({parse}: {exc})") from None
    raise JobError(f"{what}: unknown parse method {parse!r}")


def _source_text(spec: dict, run_dir: str, stdout: str) -> str:
    source = spec.get("source", "file" if "path" in spec else "stdout")
    if source == "stdout":
        return stdout
    if source == "file":
        path = os.path.join(run_dir, spec["path"])
        if not os.path.exists(path):
            raise JobError(f"the file {spec['path']!r} was not written in the run directory")
        with open(path, errors="replace") as f:
            return f.read()
    raise JobError(f"unknown source {source!r} (stdout or file)")


def _read_log(path: str, warn: Callable[[str], None] | None = None) -> list[dict]:
    """The records of runs.jsonl; an unparsable (truncated) line is skipped with a warning."""
    rows: list[dict] = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if not line.strip():
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                if warn:
                    warn(f"  (skipping unreadable line {n} of {path})")
    return rows


# ----------------------------------------------------------------- the solver
def job_fingerprint(job) -> str:
    """What a cached run depends on besides its parameters: the command, how
    the observable and the cost are read, the failure rules, the environment and the
    content of the template files.  A record whose fingerprint differs is not reused."""
    h = hashlib.sha256()
    # "solver" is always None now; it stays so that run logs written by 2.0.0 keep matching
    parts = {"command": job.command, "solver": None, "observable": job.observable,
             "cost": job.cost, "failure": job.failure, "env": job.env,
             "render_files": job.render_files}
    h.update(json.dumps(parts, sort_keys=True, default=str).encode())
    files = []
    for src in sorted(job.files or {}):
        p = src if os.path.isabs(src) else os.path.join(job.base_dir, src)
        files.append(p)
    if job.template_dir:
        for root, _, names in os.walk(job.template_dir):
            files += [os.path.join(root, n) for n in sorted(names)]
    for p in files:
        try:
            with open(p, "rb") as f:
                h.update(os.path.relpath(p, job.base_dir).encode() + f.read())
        except OSError:
            h.update(p.encode())
    return h.hexdigest()[:16]


class CommandSolver:
    """A callable solver built from a job file.  See the module docstring.

    fresh=True ignores the run log (but still appends to it); retry_failed=True runs a
    setting again when the log holds only a failed run for it (what `errbudget check` does,
    so that a fixed command is actually tried).
    """

    def __init__(self, job, log: Callable[[str], None] | None = print, fresh: bool = False,
                 retry_failed: bool = False):
        self.job = job
        self.log = log or (lambda s: None)
        self.workdir = job.workdir
        os.makedirs(self.workdir, exist_ok=True)
        self.log_path = os.path.join(self.workdir, "runs.jsonl")
        self.fingerprint = job_fingerprint(job)
        self.retry_failed = retry_failed
        self.cache: dict[str, dict] = {}
        self.n_stale = 0
        rows = _read_log(self.log_path, warn=self.log)
        self.next_id = 1 + max((int(r.get("run_id", 0)) for r in rows), default=0)
        if not fresh:
            for r in rows:
                if r.get("fingerprint") == self.fingerprint:
                    self.cache[r["key"]] = r
                else:
                    self.n_stale += 1
            if self.n_stale:
                self.log(f"  ({self.n_stale} logged runs belong to an earlier version of the "
                         f"job file and are not reused)")
        self.n_calls = 0
        self.n_cached = 0

    # -- helpers
    @staticmethod
    def key_of(params: dict) -> str:
        return json.dumps({k: params[k] for k in sorted(params)}, sort_keys=True, default=str)

    def _record(self, rec: dict) -> None:
        self.cache[rec["key"]] = rec
        with open(self.log_path, "a") as f:
            f.write(json.dumps(rec, default=str) + "\n")

    # -- the call the engine makes
    def __call__(self, **params):
        self.n_calls += 1
        key = self.key_of(params)
        rec = self.cache.get(key)
        if rec is not None and not (self.retry_failed and rec["status"] == "failed"):
            self.n_cached += 1
            self.log(f"  run {rec['run_id']:>4} (cached)  {self._pstr(params)}  -> "
                     f"{self._vstr(rec)}")
            return self._outcome(rec)
        run_id = f"{self.next_id:04d}"
        self.next_id += 1
        rec = {"run_id": run_id, "key": key, "params": params, "status": "running",
               "fingerprint": self.fingerprint, "started": time.strftime("%Y-%m-%d %H:%M:%S")}
        rec.update(self._run_command(run_id, params))
        self._record(rec)
        self.log(f"  run {run_id:>4}           {self._pstr(params)}  -> {self._vstr(rec)}")
        return self._outcome(rec)

    def _pstr(self, params: dict) -> str:
        return ", ".join(f"{k}={format_value(v)}" for k, v in params.items())

    @staticmethod
    def _vstr(rec: dict) -> str:
        if rec["status"] == "ok":
            return f"{rec['value']:.8g}  (cost {rec['cost']:.4g})"
        if rec["status"] == "unconverged":
            return f"{rec.get('value', float('nan')):.8g}  NOT CONVERGED"
        return f"FAILED: {rec.get('error', '')[:80]}"

    def _outcome(self, rec: dict):
        name = self.job.observable["name"]
        if rec["status"] == "ok":
            return {name: rec["value"], "converged": True}, float(rec["cost"])
        if rec["status"] == "unconverged":
            return {name: rec.get("value", float("nan")), "converged": False}, float(rec["cost"])
        raise RunFailed(rec.get("error", "failed"))

    # -- shell command
    def _run_command(self, run_id: str, params: dict) -> dict:
        job = self.job
        run_dir = os.path.abspath(os.path.join(self.workdir, f"run_{run_id}"))
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)
        try:
            command = prepare_run_dir(job, run_dir, params, run_id)
        except JobError as exc:
            return {"status": "failed", "error": str(exc), "cost": 0.0, "run_dir": run_dir}
        env = dict(os.environ, **{k: str(v) for k, v in (job.env or {}).items()})
        t0 = time.perf_counter()
        code, timed_out = None, False
        proc = subprocess.Popen(command, shell=True, cwd=run_dir, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                start_new_session=True)
        try:
            out_b, err_b = proc.communicate(timeout=job.timeout)
            code = proc.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
            try:                                     # the whole process group, not only /bin/sh
                os.killpg(proc.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass
            out_b, err_b = proc.communicate()
        stdout = out_b.decode(errors="replace") if out_b else ""
        stderr = err_b.decode(errors="replace") if err_b else ""
        wall = time.perf_counter() - t0
        with open(os.path.join(run_dir, "stdout.txt"), "w") as f:
            f.write(stdout)
        with open(os.path.join(run_dir, "stderr.txt"), "w") as f:
            f.write(stderr)
        rec = {"run_dir": run_dir, "exit_code": code, "wall": wall, "cost": wall,
               "command": command}
        fail = job.failure or {}
        if timed_out:
            return dict(rec, status="failed", error=f"timed out after {job.timeout} s")
        if code != 0 and not fail.get("allow_nonzero_exit", False):
            tail = (stderr.strip().splitlines() or stdout.strip().splitlines() or [""])[-1]
            return dict(rec, status="failed", error=f"exit code {code}: {tail[:200]}")
        # the texts the failure patterns are searched in: stdout, stderr, the files the
        # observable and the cost are read from, and failure.path if given
        texts = [stdout, stderr]
        for spec in (job.observable, job.cost):
            if spec and "path" in spec:
                try:
                    texts.append(_source_text(spec, run_dir, stdout))
                except JobError:
                    pass
        for p in _as_list(fail.get("path")):
            fp = os.path.join(run_dir, p)
            if os.path.exists(fp):
                with open(fp, errors="replace") as f:
                    texts.append(f.read())

        def found(pat: str) -> bool:
            return any(re.search(pat, t, re.MULTILINE) for t in texts)

        for pat in _as_list(fail.get("failed_if_regex")):
            if found(pat):
                return dict(rec, status="failed", error=f"output matched {pat!r}")
        try:
            value = extract_number(job.observable, _source_text(job.observable, run_dir, stdout),
                                   f"observable {job.observable['name']!r}")
        except JobError as exc:
            with open(os.path.join(run_dir, "error.txt"), "w") as f:
                f.write(str(exc) + "\n")
            return dict(rec, status="failed", error=str(exc))
        if job.cost:
            try:
                rec["cost"] = extract_number(job.cost, _source_text(job.cost, run_dir, stdout), "cost")
            except JobError as exc:
                return dict(rec, status="failed", error=str(exc))
        status = "ok"
        for pat in _as_list(fail.get("unconverged_if_regex")):
            if found(pat):
                status = "unconverged"
        pat = fail.get("converged_if_regex")
        if pat and not found(pat):
            status = "unconverged"
        if not job.keep_runs:
            shutil.rmtree(run_dir, ignore_errors=True)
        return dict(rec, status=status, value=value)


def prepare_run_dir(job, run_dir: str, params: dict, run_id: str = "0000") -> str:
    """Create the run directory (template copied, files rendered, params.json written)
    and return the rendered command.  Raises JobError for an unknown placeholder."""
    if job.template_dir:
        shutil.copytree(job.template_dir, run_dir)
    else:
        os.makedirs(run_dir, exist_ok=True)
    values = dict(params, run_dir=run_dir, run_id=run_id, job_dir=job.base_dir)
    with open(os.path.join(run_dir, "params.json"), "w") as f:
        json.dump(params, f, indent=1, default=str)
    targets = []
    for src, dst in (job.files or {}).items():
        src_path = src if os.path.isabs(src) else os.path.join(job.base_dir, src)
        dst_path = os.path.join(run_dir, dst)
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        shutil.copy(src_path, dst_path)
        targets.append(dst_path)
    if job.render_files is not None:
        targets += [os.path.join(run_dir, p) for p in job.render_files]
    elif job.template_dir:
        for root, _, files in os.walk(run_dir):
            for fn in files:
                p = os.path.join(root, fn)
                if os.path.getsize(p) < 4_000_000:
                    targets.append(p)
    for p in dict.fromkeys(targets):
        try:
            with open(p, encoding="utf-8") as f:
                text = f.read()
        except (UnicodeDecodeError, OSError):
            continue
        if "{{" in text:
            with open(p, "w", encoding="utf-8") as f:
                f.write(render(text, values, os.path.relpath(p, run_dir)))
    return render(job.command, values, "command") if job.command else ""


def dry_render(job, params: dict) -> str:
    """Render the command and the template files for `params` in a temporary directory
    (nothing is run); returns the rendered command.  Used by `errbudget check`."""
    tmp = tempfile.mkdtemp(prefix="errbudget_check_")
    try:
        return prepare_run_dir(job, os.path.join(tmp, "run"), params, "check")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
