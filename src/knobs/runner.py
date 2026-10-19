"""Run a solver that is a shell command (or an importable Python function) from a job file.

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

import importlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time
from typing import Any, Callable

from .api import RunFailed

PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\s*(?::([^}]*))?\s*\}\}")
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")


class JobError(ValueError):
    """A job file that cannot be used as written; the message says what to fix."""


def _as_list(x) -> list:
    return [] if x is None else [x] if isinstance(x, str) else list(x)


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
        return format_value(values[name], fmt)
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


def extract_number(spec: dict, text: str, what: str = "observable") -> float:
    """Read one number out of `text` as the spec says.

    spec["parse"]:
      "regex"        spec["regex"] with one capture group; the last match (or the
                     first with spec["which"] = "first")
      "json"         json.loads(text) and spec["key"] as a path ("forces.drag", "a[0].b")
      "column"       the last non-comment line split on whitespace or commas,
                     column spec["column"] (1-based); spec["comment"] holds the
                     characters that start a comment line (default '#'; '#@' for xvg)
      "last_number"  the last number in the text (default when nothing else is given)
    """
    parse = spec.get("parse")
    if parse is None:
        parse = "regex" if "regex" in spec else "json" if "key" in spec else \
            "column" if "column" in spec else "last_number"
    try:
        if parse == "regex":
            ms = re.findall(spec["regex"], text, re.MULTILINE)
            if not ms:
                raise JobError(f"{what}: the pattern {spec['regex']!r} matched nothing")
            m = ms[0] if spec.get("which") == "first" else ms[-1]
            if isinstance(m, tuple):
                m = m[0]
            return float(str(m).replace("D", "e").replace("d", "e"))
        if parse == "json":
            try:
                obj = json.loads(text)
            except json.JSONDecodeError:
                # the output has other lines around the JSON: take the last line that parses
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
            return float(_json_path(obj, spec["key"]))
        if parse == "column":
            comment = str(spec.get("comment", "#"))          # characters that start a comment line
            lines = [ln for ln in text.splitlines() if ln.strip() and ln.lstrip()[0] not in comment]
            if not lines:
                raise JobError(f"{what}: no data lines")
            row = re.split(r"[,\s]+", lines[-1].strip())
            return float(row[int(spec["column"]) - 1])
        if parse == "last_number":
            ms = NUMBER.findall(text)
            if not ms:
                raise JobError(f"{what}: no number found in the output")
            return float(ms[-1].replace("D", "e").replace("d", "e"))
    except (KeyError, IndexError, ValueError, json.JSONDecodeError) as exc:
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


# ----------------------------------------------------------------- the solver
def load_python_solver(ref: str, base_dir: str) -> Callable:
    """'package.module:function' or 'path/to/file.py:function'."""
    if ":" not in ref:
        raise JobError(f"solver must be 'module:function' or 'file.py:function', got {ref!r}")
    mod, func = ref.rsplit(":", 1)
    if mod.endswith(".py"):
        path = mod if os.path.isabs(mod) else os.path.join(base_dir, mod)
        spec = importlib.util.spec_from_file_location("knobs_user_solver", path)
        module = importlib.util.module_from_spec(spec)
        sys.modules["knobs_user_solver"] = module
        spec.loader.exec_module(module)
    else:
        if base_dir not in sys.path:
            sys.path.insert(0, base_dir)
        module = importlib.import_module(mod)
    try:
        return getattr(module, func)
    except AttributeError:
        raise JobError(f"{mod} has no function {func!r}") from None


class CommandSolver:
    """A callable solver built from a job file.  See the module docstring."""

    def __init__(self, job, log: Callable[[str], None] | None = print, fresh: bool = False):
        self.job = job
        self.log = log or (lambda s: None)
        self.workdir = job.workdir
        os.makedirs(self.workdir, exist_ok=True)
        self.log_path = os.path.join(self.workdir, "runs.jsonl")
        self.cache: dict[str, dict] = {}
        self.n_done = 0
        if not fresh and os.path.exists(self.log_path):
            with open(self.log_path) as f:
                for line in f:
                    if line.strip():
                        r = json.loads(line)
                        self.cache[r["key"]] = r
            self.n_done = len(self.cache)
        self.n_calls = 0
        self.n_cached = 0
        self._py = load_python_solver(job.solver, job.base_dir) if job.solver else None

    # -- helpers
    @staticmethod
    def key_of(params: dict) -> str:
        return json.dumps({k: params[k] for k in sorted(params)}, sort_keys=True, default=str)

    def _record(self, rec: dict) -> None:
        self.cache[rec["key"]] = rec
        with open(self.log_path, "a") as f:
            f.write(json.dumps(rec, default=str) + "\n")

    def _say(self, msg: str) -> None:
        self.log(msg)

    # -- the call the engine makes
    def __call__(self, **params):
        self.n_calls += 1
        key = self.key_of(params)
        if key in self.cache:
            rec = self.cache[key]
            self.n_cached += 1
            self._say(f"  run {rec['run_id']:>4} (cached)  {self._pstr(params)}  -> "
                      f"{self._vstr(rec)}")
            return self._outcome(rec)
        run_id = f"{len(self.cache) + 1:04d}"
        rec = {"run_id": run_id, "key": key, "params": params, "status": "running",
               "started": time.strftime("%Y-%m-%d %H:%M:%S")}
        if self._py is not None:
            rec.update(self._run_python(params))
        else:
            rec.update(self._run_command(run_id, params))
        self._record(rec)
        self._say(f"  run {run_id:>4}           {self._pstr(params)}  -> {self._vstr(rec)}")
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

    # -- python function
    def _run_python(self, params: dict) -> dict:
        t0 = time.perf_counter()
        try:
            out = self._py(**params)
        except RunFailed as exc:
            return {"status": "failed", "error": str(exc) or "RunFailed",
                    "cost": time.perf_counter() - t0}
        wall = time.perf_counter() - t0
        cost = wall
        if isinstance(out, tuple) and len(out) == 2:
            out, cost = out
        converged = True
        if isinstance(out, dict):
            out = dict(out)
            converged = bool(out.pop("converged", True))
            name = self.job.observable["name"]
            if name in out:
                value = float(out[name])
            elif len(out) == 1:
                value = float(next(iter(out.values())))
            else:
                return {"status": "failed", "cost": float(cost),
                        "error": f"solver returned {sorted(out)}, none named {name!r}"}
        else:
            value = float(out)
        return {"status": "ok" if converged else "unconverged", "value": value,
                "cost": float(cost), "wall": wall}

    # -- shell command
    def _run_command(self, run_id: str, params: dict) -> dict:
        job = self.job
        run_dir = os.path.abspath(os.path.join(self.workdir, f"run_{run_id}"))
        if os.path.exists(run_dir):
            shutil.rmtree(run_dir)
        if job.template_dir:
            shutil.copytree(job.template_dir, run_dir)
        else:
            os.makedirs(run_dir)
        values = dict(params, run_dir=run_dir, run_id=run_id, job_dir=job.base_dir)
        with open(os.path.join(run_dir, "params.json"), "w") as f:
            json.dump(params, f, indent=1, default=str)
        try:
            self._render_files(run_dir, values)
            command = render(job.command, values, "command")
        except JobError as exc:
            return {"status": "failed", "error": str(exc), "cost": 0.0, "run_dir": run_dir}
        env = dict(os.environ, **{k: str(v) for k, v in (job.env or {}).items()})
        t0 = time.perf_counter()
        try:
            proc = subprocess.run(command, shell=True, cwd=run_dir, env=env,
                                  capture_output=True, text=True, errors="replace",
                                  timeout=job.timeout)
            stdout, stderr, code, timed_out = proc.stdout, proc.stderr, proc.returncode, False
        except subprocess.TimeoutExpired as exc:
            stdout = (exc.stdout or b"").decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            stderr = (exc.stderr or b"").decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
            code, timed_out = None, True
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
        for pat in _as_list(fail.get("failed_if_regex")):
            if re.search(pat, stdout, re.MULTILINE) or re.search(pat, stderr, re.MULTILINE):
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
            if re.search(pat, stdout, re.MULTILINE) or re.search(pat, stderr, re.MULTILINE):
                status = "unconverged"
        pat = fail.get("converged_if_regex")
        if pat and not (re.search(pat, stdout, re.MULTILINE) or re.search(pat, stderr, re.MULTILINE)):
            status = "unconverged"
        if not job.keep_runs:
            shutil.rmtree(run_dir, ignore_errors=True)
        return dict(rec, status=status, value=value)

    def _render_files(self, run_dir: str, values: dict) -> None:
        job = self.job
        targets = []
        for src, dst in (job.files or {}).items():
            src_path = src if os.path.isabs(src) else os.path.join(job.base_dir, src)
            dst_path = os.path.join(run_dir, dst)
            os.makedirs(os.path.dirname(dst_path) or run_dir, exist_ok=True)
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
