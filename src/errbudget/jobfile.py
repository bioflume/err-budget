"""Job files: a tuning described in YAML, JSON or TOML instead of Python.

    errbudget init             # writes an annotated job file to start from
    errbudget check job.yaml   # validates it and runs your starting setting once
    errbudget run job.yaml     # tunes, prints the report, writes the result files

A job file says what to run (`command`, or `solver` for a Python function), which
parameters to tune (`parameters`, each with its kind and candidate values), what to
read from the output (`observable`) and how accurate it has to be (`target`).  The
reference is in docs/job-files.md; `errbudget init` writes one with every field explained.
"""
from __future__ import annotations

import inspect
import json
import os
import re
import time
from dataclasses import dataclass, field
from typing import Any, Callable

from .api import Parameter, TuneResult, tune, KINDS, CONFIGURATIONS
from ._core import run_protocol
from .runner import (BUILTIN_PLACEHOLDERS, CommandSolver, JobError, placeholders, format_value, dry_render,
                     parse_method)

KIND_ALIASES = {
    "discretisation": "discretisation", "discretization": "discretisation",
    "mesh": "discretisation", "resolution": "discretisation", "d": "discretisation",
    "nested-solver": "nested-solver", "nested_solver": "nested-solver", "nested": "nested-solver",
    "tolerance": "nested-solver", "solver": "nested-solver", "n": "nested-solver",
    "cost-only": "cost-only", "cost_only": "cost-only", "cost": "cost-only", "c": "cost-only",
}
TUNE_OPTIONS = ("screen", "accept_unconverged", "floor_probe", "configuration")
FAILURE_RULES = ("failed_if_regex", "unconverged_if_regex", "converged_if_regex",
                 "allow_nonzero_exit", "path")


_NUM = re.compile(r"^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$")


def _num(v):
    """YAML 1.1 reads 1e-4 as a string (no dot): turn numeric strings into numbers."""
    if isinstance(v, str) and _NUM.match(v.strip()):
        s = v.strip()
        if re.match(r"^[-+]?\d+$", s):
            return int(s)
        return float(s)
    return v


@dataclass
class JobSpec:
    """A validated job file.  Build one with `load_job` or `JobSpec.from_dict`."""
    parameters: list                      # [{"name", "kind", "values", "start", "h", "model"}]
    target: float
    observable: dict                      # {"name", "source", "path", "parse", "regex"/"key"/"column", ...}
    command: str | None = None
    solver: str | None = None             # "module:function" or "file.py:function"
    name: str = "errbudget job"
    relative: bool = False
    start: dict = field(default_factory=dict)
    options: dict = field(default_factory=dict)
    workdir: str = "errbudget_runs"
    template_dir: str | None = None
    files: dict = field(default_factory=dict)
    render_files: list | None = None
    env: dict = field(default_factory=dict)
    timeout: float | None = None
    keep_runs: bool = True
    cost: dict | None = None
    failure: dict = field(default_factory=dict)
    output: str = "errbudget_result"
    base_dir: str = "."
    source_path: str | None = None

    # ---------------------------------------------------------------- loading
    @classmethod
    def from_dict(cls, d: dict, base_dir: str = ".", source_path: str | None = None) -> "JobSpec":
        errors: list[str] = []
        d = dict(d)
        base_dir = os.path.abspath(base_dir)
        for key in ("target", "timeout"):
            if key in d:
                d[key] = _num(d[key])
        if isinstance(d.get("start"), dict):
            d["start"] = {k: _num(v) for k, v in d["start"].items()}
        if isinstance(d.get("parameters"), list):
            d["parameters"] = [dict(p, values=[_num(v) for v in p["values"]],
                                    **({"start": _num(p["start"])} if "start" in p else {}))
                               if isinstance(p, dict) and isinstance(p.get("values"), list) else p
                               for p in d["parameters"]]
        if isinstance(d.get("options"), dict) and isinstance(d["options"].get("floor_probe"), list):
            d["options"] = dict(d["options"], floor_probe=[
                {k: _num(v) for k, v in x.items()} if isinstance(x, dict) else x
                for x in d["options"]["floor_probe"]])

        def need(key, typ):
            if key not in d:
                errors.append(f"missing field {key!r}")
                return None
            if not isinstance(d[key], typ):
                errors.append(f"{key!r} must be {typ.__name__ if not isinstance(typ, tuple) else '/'.join(t.__name__ for t in typ)}")
                return None
            return d[key]

        command, solver = d.get("command"), d.get("solver")
        if isinstance(command, list):
            command = " ".join(str(c) for c in command)
        if not command and not solver:
            errors.append("give either 'command' (a shell command with {{placeholders}}) or "
                          "'solver' (module:function)")
        if command and solver:
            errors.append("give 'command' or 'solver', not both")

        params_in = need("parameters", list) or []
        params, names = [], []
        for i, p in enumerate(params_in):
            if not isinstance(p, dict) or "name" not in p:
                errors.append(f"parameters[{i}] needs at least 'name' and 'values'")
                continue
            kind = KIND_ALIASES.get(str(p.get("kind", "discretisation")).lower())
            if kind is None:
                errors.append(f"parameter {p['name']!r}: kind must be one of {KINDS} "
                              f"(got {p.get('kind')!r})")
            vals = p.get("values")
            if not isinstance(vals, list) or len(vals) < 2:
                errors.append(f"parameter {p['name']!r}: 'values' must list at least two "
                              f"candidate values, coarse (cheap) to fine (accurate)")
                vals = vals if isinstance(vals, list) else []
            q = dict(name=str(p["name"]), kind=kind, values=vals,
                     model=p.get("model", "power"), h=p.get("h"), start=p.get("start"))
            if q["model"] not in ("power", "exponential"):
                errors.append(f"parameter {p['name']!r}: model must be 'power' or 'exponential'")
            hh = q["h"]
            if hh is not None and not (hh in ("value", "inverse") or
                                       (isinstance(hh, list) and len(hh) == len(vals)
                                        and all(isinstance(_num(x), (int, float)) for x in hh))):
                errors.append(f"parameter {p['name']!r}: h must be 'value', 'inverse' or a list "
                              f"of {len(vals)} numbers (one per value)")
            elif isinstance(hh, list):
                q["h"] = [float(_num(x)) for x in hh]
            if q["start"] is not None and q["start"] not in vals:
                errors.append(f"parameter {p['name']!r}: start value {q['start']!r} is not "
                              f"among its values {vals!r}")
            if q["name"] in names:
                errors.append(f"parameter {q['name']!r} is listed twice")
            names.append(q["name"])
            params.append(q)
        if not params and not errors:
            errors.append("'parameters' is empty")

        target = need("target", (int, float))
        if target is not None and target <= 0:
            errors.append("'target' must be positive")

        obs = d.get("observable")
        if isinstance(obs, str):
            obs = {"name": obs}
        if not isinstance(obs, dict) or "name" not in obs:
            errors.append("'observable' needs at least a 'name' (and how to read it: "
                          "regex, key, column, or the last number of stdout)")
            obs = obs if isinstance(obs, dict) else {}
        if obs.get("source") == "file" and "path" not in obs:
            errors.append("observable: 'source: file' needs 'path'")

        start = dict(d.get("start") or {})
        for p in params:
            if p["start"] is not None:
                start.setdefault(p["name"], p["start"])
        for k, v in start.items():
            p = next((q for q in params if q["name"] == k), None)
            if p is None:
                errors.append(f"start mentions an unknown parameter {k!r}")
            elif v not in p["values"]:
                errors.append(f"start value {v!r} for {k!r} is not among its values {p['values']!r}")

        options = dict(d.get("options") or {})
        cfg = str(options.get("configuration", "v1.2.1"))
        if cfg not in CONFIGURATIONS:
            errors.append("options.configuration must be v1.2.1 (default), v1.2 or v1.1")
        for key in ("screen", "accept_unconverged"):
            if key in options and not isinstance(options[key], bool):
                errors.append(f"options.{key} must be true or false")
        fp = options.get("floor_probe")
        if fp is not None and (not isinstance(fp, list) or len(fp) < 2
                               or not all(isinstance(x, dict) for x in fp)):
            errors.append("options.floor_probe must list two or more settings (dicts)")
        engine_options = set(inspect.signature(run_protocol).parameters) - {
            "lookup", "ladders", "anchor", "eps", "hs", "kinds", "models", "noise_floor", "use_screen"}
        for key in list(options):
            if key in TUNE_OPTIONS:
                continue
            if key not in engine_options:
                errors.append(f"options.{key} is not an option; the engine options are "
                              f"{', '.join(sorted(engine_options))}")
            else:
                options[key] = _num(options[key])

        timeout = d.get("timeout")
        if timeout is not None and not (isinstance(timeout, (int, float)) and not isinstance(timeout, bool)
                                        and timeout > 0):
            errors.append("'timeout' must be a positive number of seconds")

        if command:
            known = set(names) | set(BUILTIN_PLACEHOLDERS)
            for ph in placeholders(command):
                if ph not in known:
                    errors.append(f"command uses {{{{{ph}}}}}, which is not a parameter")

        tdir = d.get("template_dir")
        if tdir is not None:
            tdir = tdir if os.path.isabs(tdir) else os.path.join(base_dir, tdir)
            if not os.path.isdir(tdir):
                errors.append(f"template_dir {d['template_dir']!r} is not a directory")
        files = d.get("files") or {}
        if not isinstance(files, dict):
            errors.append("'files' must map source file -> destination in the run directory")
            files = {}
        for src in files:
            p = src if os.path.isabs(src) else os.path.join(base_dir, src)
            if not os.path.isfile(p):
                errors.append(f"files: {src!r} does not exist")

        cost = d.get("cost")
        if isinstance(cost, str):
            cost = None if cost in ("wall", "time", "wall_time") else {"regex": cost}
        if cost is not None and not isinstance(cost, dict):
            errors.append("'cost' must be 'wall' (default) or a dict saying how to read it")
            cost = None
        if cost is not None and cost.get("source") == "file" and "path" not in cost:
            errors.append("cost: 'source: file' needs 'path'")
        failure = dict(d.get("failure") or {})
        if not isinstance(d.get("failure") or {}, dict):
            errors.append(f"'failure' must be a mapping ({', '.join(FAILURE_RULES)})")
            failure = {}
        for key in set(failure) - set(FAILURE_RULES):
            errors.append(f"failure.{key} is not a failure rule")
        # every regular expression must compile
        for what, spec in (("observable", obs), ("cost", cost or {})):
            if "regex" in spec:
                try:
                    if re.compile(str(spec["regex"])).groups < 1:
                        errors.append(f"{what}.regex needs one capture group around the number")
                except re.error as exc:
                    errors.append(f"{what}.regex does not compile: {exc}")
            if parse_method(spec) not in ("regex", "json", "column", "last_number"):
                errors.append(f"{what}.parse must be regex, json, column or last_number")
        for key in ("failed_if_regex", "unconverged_if_regex", "converged_if_regex"):
            pats = failure.get(key)
            for pat in ([pats] if isinstance(pats, str) else pats or []):
                try:
                    re.compile(str(pat))
                except re.error as exc:
                    errors.append(f"failure.{key} {pat!r} does not compile: {exc}")
        if solver and ":" in str(solver):
            mod = str(solver).rsplit(":", 1)[0]
            if mod.endswith(".py") and not os.path.isfile(mod if os.path.isabs(mod)
                                                          else os.path.join(base_dir, mod)):
                errors.append(f"solver file {mod!r} does not exist")
        elif solver:
            errors.append("solver must be 'module:function' or 'file.py:function'")
        unknown = set(d) - {"name", "command", "solver", "parameters", "target", "observable",
                            "relative", "start", "options", "workdir", "template_dir", "files",
                            "render_files", "env", "timeout", "keep_runs", "cost", "failure",
                            "output"}
        if unknown:
            errors.append(f"unknown field(s): {', '.join(sorted(unknown))}")
        if errors:
            raise JobError("the job file cannot be used as written:\n  - " + "\n  - ".join(errors))

        workdir = d.get("workdir", "errbudget_runs")
        workdir = workdir if os.path.isabs(workdir) else os.path.join(base_dir, workdir)
        output = d.get("output", "errbudget_result")
        output = output if os.path.isabs(output) else os.path.join(base_dir, output)
        job = cls(parameters=params, target=float(target), observable=obs, command=command,
                  solver=solver, name=str(d.get("name", "errbudget job")),
                  relative=bool(d.get("relative", False)), start=start, options=options,
                  workdir=workdir, template_dir=tdir, files=files,
                  render_files=d.get("render_files"), env=dict(d.get("env") or {}),
                  timeout=timeout, keep_runs=bool(d.get("keep_runs", True)),
                  cost=cost, failure=failure, output=output, base_dir=base_dir,
                  source_path=source_path)
        try:                       # the engine's own checks on the candidate lists (direction, h)
            job.to_parameters()
        except (ValueError, TypeError, KeyError) as exc:
            raise JobError(f"the job file cannot be used as written:\n  - {exc}") from None
        return job

    # ---------------------------------------------------------------- building
    def to_parameters(self) -> list[Parameter]:
        out = []
        for p in self.parameters:
            h = None
            if p.get("h") == "value":
                h = lambda v: float(v)
            elif p.get("h") == "inverse":
                h = lambda v: 1.0 / float(v)
            elif isinstance(p.get("h"), list):
                table = {v: float(x) for v, x in zip(p["values"], p["h"])}
                h = lambda v, table=table: table[v]
            out.append(Parameter(p["name"], p["values"], kind=p["kind"], h=h, model=p.get("model", "power")))
        return out

    def tune_kwargs(self) -> dict:
        o = self.options
        kw: dict[str, Any] = dict(target=self.target, observable=self.observable["name"],
                                  start=self.start or None, relative=self.relative,
                                  screen=bool(o.get("screen", False)),
                                  accept_unconverged=bool(o.get("accept_unconverged", False)),
                                  floor_probe=o.get("floor_probe"))
        kw.update(CONFIGURATIONS[str(o.get("configuration", "v1.2.1"))])
        for key, val in o.items():
            if key not in TUNE_OPTIONS:
                kw[key] = val           # expert overrides, passed to the engine as given
        return kw

    def describe(self) -> str:
        """A short text description: what runs, the parameters, the target."""
        L = [f"job: {self.name}"]
        L.append(f"runs: {self.command!r}" if self.command else f"solver: {self.solver}")
        L.append(f"observable: {self.observable['name']} "
                 f"({self.observable.get('parse') or parse_method(self.observable).replace('_', ' ')}"
                 f" from {self.observable.get('source', 'file' if 'path' in self.observable else 'stdout')})")
        L.append(f"target: {self.target:g}{' relative' if self.relative else ' absolute'}")
        L.append("parameters (coarse -> fine):")
        for p in self.parameters:
            vals = ", ".join(format_value(v) for v in p["values"])
            st = self.start.get(p["name"], p["values"][len(p["values"]) // 2])
            L.append(f"  {p['name']:<18} {p['kind']:<15} [{vals}]   start {format_value(st)}")
        L.append(f"configuration: {self.options.get('configuration', 'v1.2.1')}"
                 f"{'; preliminary screen' if self.options.get('screen') else ''}"
                 f"{'; unconverged runs accepted' if self.options.get('accept_unconverged') else ''}")
        L.append(f"run directories: {self.workdir}")
        return "\n".join(L)


def load_job(path: str) -> JobSpec:
    """Read a YAML (.yaml/.yml), JSON (.json) or TOML (.toml) job file."""
    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise JobError(f"no such job file: {path}")
    ext = os.path.splitext(path)[1].lower()
    with open(path, "rb") as f:
        raw = f.read()
    if ext in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError:
            raise JobError("reading YAML needs PyYAML: pip install pyyaml  (or write the job "
                           "file as JSON or TOML)") from None
        d = yaml.safe_load(raw)
    elif ext == ".json":
        d = json.loads(raw)
    elif ext == ".toml":
        try:
            import tomllib
        except ImportError:                                   # Python < 3.11
            raise JobError("reading TOML needs Python 3.11+ (or write the job file as JSON)") from None
        d = tomllib.loads(raw.decode())
    else:
        raise JobError(f"unknown job file type {ext!r}: use .yaml, .json or .toml")
    if not isinstance(d, dict):
        raise JobError("the job file must be a mapping of fields")
    return JobSpec.from_dict(d, base_dir=os.path.dirname(path), source_path=path)


# ------------------------------------------------------------------- running
def run_job(job: JobSpec, log: Callable[[str], None] | None = print, fresh: bool = False,
            write: bool = True) -> TuneResult:
    """Tune the job's solver; write <output>.json and <output>.md unless write=False."""
    solver = CommandSolver(job, log=log, fresh=fresh)
    params = job.to_parameters()
    t0 = time.time()
    result = tune(solver, params, **job.tune_kwargs())
    elapsed = time.time() - t0
    if write:
        os.makedirs(os.path.dirname(os.path.abspath(job.output)), exist_ok=True)
        result.to_json(job.output + ".json")
        with open(job.output + ".md", "w") as f:
            f.write(f"# {job.name}\n\n")
            f.write(f"Job file: `{job.source_path or '-'}`; {solver.n_calls} solver calls, "
                    f"{solver.n_cached} answered from the run log; {elapsed:.0f} s.\n\n")
            f.write(result.markdown())
    return result


def check_job(job: JobSpec, log: Callable[[str], None] = print, run: bool = True) -> bool:
    """Print the job and run the starting setting once.  Returns True if that run gave a number."""
    log(job.describe())
    start = {p["name"]: job.start.get(p["name"], p["values"][len(p["values"]) // 2])
             for p in job.parameters}
    if job.command:
        cmd = dry_render(job, start)           # unknown placeholders in templates show up here
        log(f"command at the starting setting: {cmd}")
    if not run:
        log("(validation only; nothing was run)")
        return True
    log("")
    log("running the starting setting once ...")
    solver = CommandSolver(job, log=log, fresh=False, retry_failed=True)
    try:
        out, cost = solver(**start)
    except Exception as exc:                                   # noqa: BLE001
        log(f"FAILED: {exc}")
        rec = solver.cache.get(solver.key_of(start), {})
        if rec.get("run_dir"):
            log(f"look in {rec['run_dir']} (stdout.txt, stderr.txt, params.json)")
        return False
    name = job.observable["name"]
    log(f"{name} = {out[name]:.10g}   cost = {cost:.6g}"
        f"{'' if out.get('converged', True) else '   (reported as not converged)'}")
    if job.relative:
        log(f"relative target {job.target:g} means an absolute error of {job.target * abs(out[name]):.3g}")
    log("the job file works; next:  errbudget run " + (job.source_path or "job.yaml"))
    return True


TEMPLATE_YAML = '''# errbudget job file -- edit the five numbered parts, then:  errbudget check {fname}  and  errbudget run {fname}
name: my simulation

# 1. The command that runs ONE simulation.  {{{{name}}}} is replaced by a parameter value.
#    The command runs inside a fresh directory made for that run ({{{{run_dir}}}}), so refer
#    to scripts and inputs next to this file as {{{{job_dir}}}}/...
#    Examples:
#      python3 {{{{job_dir}}}}/solver.py --steps {{{{steps}}}} --tol {{{{tol}}}}
#      lmp -in in.lammps -var cutoff {{{{cutoff}}}} -var order {{{{order}}}} > log.lammps
#      blockMesh > log.mesh && simpleFoam > log.simpleFoam && tail -1 postProcessing/forceCoeffs/0/coefficient.dat
command: python3 {{{{job_dir}}}}/solver.py --steps {{{{steps}}}} --tol {{{{tol}}}} --refresh {{{{refresh}}}}

# Files copied into every run directory with their placeholders filled in.  Use
# template_dir for a whole case (an OpenFOAM case, a LAMMPS input deck); every text file
# in it that contains {{{{...}}}} is rendered.
# template_dir: case
# files: {{"templates/in.lammps": "in.lammps"}}

# 2. The parameters to tune.  Values run COARSE (cheap, inaccurate) to FINE.
#    kind: discretisation  -- refining improves the answer and costs more (mesh, time
#                             step, cutoff, k-points, basis size)
#          nested-solver   -- a tolerance or iteration cap of a solver inside the outer
#                             problem (Newton, linear solver, SCF, constraint tolerance)
#          cost-only       -- changes cost but not the answer (preconditioner fill,
#                             how often a Jacobian or neighbour list is rebuilt)
#    start: the setting you would otherwise have used (where the measurement begins).
parameters:
  - {{name: steps,   kind: discretisation, values: [25, 50, 100, 200, 400, 800], start: 100}}
  - {{name: tol,     kind: nested-solver,  values: [1e-2, 1e-4, 1e-6, 1e-8, 1e-10], start: 1e-10}}
  - {{name: refresh, kind: cost-only,      values: [16, 8, 4, 2, 1], start: 1}}
#   model: exponential   -- add to a parameter whose error decays exponentially (a
#                           plane-wave or spectral cutoff; its values must increase)

# 3. What to read from the output, and from where.
observable:
  name: u_final
  source: stdout              # stdout, or file (then give path: relative to the run directory)
  regex: "u_final\\\\s*=\\\\s*([-+0-9.eE]+)"   # first capture group; the last match is used
  # key: forces.drag          # instead of regex: a JSON output and a key path
  # column: 2                 # instead: the last data line of a table, column 2 (1-based)
  # (with none of these, the last number in the text is used)

# 4. The accuracy you need in that observable: absolute, or relative to its value.
target: 1e-4
relative: false

# How the cost of a run is measured: wall-clock time (default), or a number read from
# the output like the observable (e.g. OpenFOAM's "ExecutionTime = 12.3 s").
# cost: {{regex: "ExecutionTime\\\\s*=\\\\s*([0-9.]+)"}}
# cost: {{regex: "work\\\\s*=\\\\s*([0-9.eE+-]+)"}}

# When is a run a failure?  A non-zero exit code is (unless allow_nonzero_exit: true).
# Add patterns that mark a crash (failed_if_regex) or a run whose inner solver did not
# converge (unconverged_if_regex, or converged_if_regex for a line that must appear).
# Patterns are searched in stdout, stderr, the files the observable and the cost are read
# from, and any file named in path.
failure:
  failed_if_regex: ["FOAM FATAL", "Traceback"]
  # unconverged_if_regex: "did not converge"
  # converged_if_regex: "SIMPLE solution converged"
  # path: log.simpleFoam

# 5. Options.  configuration v1.2.1 is the article's default; screen runs a cheap
#    preliminary sweep (fewer runs, more expensive recommendations).
options:
  configuration: v1.2.1
  screen: false
  accept_unconverged: false
  # floor_probe:               # two settings refined beyond every value: measures the
  #   - {{steps: 1600, tol: 1e-12, refresh: 1}}    # round-off floor; targets near it are
  #   - {{steps: 3200, tol: 1e-12, refresh: 1}}    # declined instead of searched for

workdir: errbudget_runs        # one directory per run, plus runs.jsonl (the log; a stopped
                               # tuning resumes from it)
output: errbudget_result       # writes errbudget_result.json and errbudget_result.md
timeout: 3600                  # seconds per run
keep_runs: true
'''


def init_job(path: str = "errbudget.yaml", overwrite: bool = False) -> str:
    """Write an annotated job file; returns the path."""
    if os.path.exists(path) and not overwrite:
        raise JobError(f"{path} exists; choose another name or pass --force")
    with open(path, "w") as f:
        f.write(TEMPLATE_YAML.format(fname=os.path.basename(path)))
    return path
