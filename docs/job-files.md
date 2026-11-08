# Job files: the reference

A job file describes one tuning: what to run, which parameters to vary, what to read
from the output, and how accurate it has to be. It is YAML (`.yaml`, needs PyYAML:
`pip install "errbudget[yaml]"`), JSON (`.json`) or TOML (`.toml`, Python 3.11+). `errbudget init`
writes an annotated one to start from.

```bash
errbudget init job.yaml   # template
errbudget check job.yaml  # validate; run the starting setting once and show what was read
errbudget run job.yaml    # tune
errbudget runs job.yaml   # the run log
errbudget report errbudget_result.json
```

## The fields

### `command` (required)

`command` is a shell command that runs **one** simulation. `{{name}}` is replaced by the
value of the parameter `name`; three more placeholders are always available:

| placeholder | value |
|---|---|
| `{{run_dir}}` | the directory made for this run, in which the command runs |
| `{{run_id}}` | the run's number (`0001`, `0002`, …) |
| `{{job_dir}}` | the directory of the job file |

The command runs **inside the run directory**, so a script or input next to the job file
is `{{job_dir}}/solver.py`. A format spec can follow the name: `{{dt:.3e}}`. Integers are
written as integers, floats in their shortest form (`1e-06`, `0.05`). Braces that are not
`{{name}}` are left alone, so OpenFOAM dictionaries are safe.

### `parameters` (required)

A list; each entry has

| key | |
|---|---|
| `name` | the placeholder name |
| `kind` | `discretisation` (default), `nested-solver` or `cost-only`; see the README. Aliases: `mesh`, `resolution`, `tolerance`, `nested`, `cost` |
| `values` | at least two candidate values, **coarse (cheap, inaccurate) first, fine last** |
| `start` | the value you would otherwise have used; the measurement begins there (default: the middle value) |
| `model` | `power` (default) or `exponential` for a cutoff whose error decays exponentially; such values must increase |
| `h` | how the refinement scalar follows from a value when the default is wrong: `value` (the value itself decreases as it refines), `inverse` (1/value: the value increases as it refines), or a list of numbers, one per value, that decreases coarse to fine. Needed for an increasing list of values below 1, such as a cutoff radius in nm |

The order matters: the tool must know which end is finer. Counts (mesh points, time steps)
increase as they refine; tolerances decrease. An increasing list of values below 1 is
refused as ambiguous: reverse it, or give `h`.

A top-level `start:` mapping may set starting values instead of the per-parameter `start`.

### `observable` (required)

What to read, and from where. Either the name alone (the last number on stdout) or

| key | |
|---|---|
| `name` | the observable's name (appears in the report) |
| `source` | `stdout` (default) or `file` |
| `path` | with `source: file`: the file, relative to the run directory |
| `regex` | a pattern with one capture group; the **last** match is used (`which: first` for the first). Write `\s` as `\\s` in a double-quoted YAML string, or use single quotes |
| `key` | the output is JSON and this is the key path: `forces.drag`, `rows[2].value`. If stdout has other lines around the JSON, the last line that is JSON is used |
| `column` | the output is a table: the last non-comment line, this column (1-based), split on whitespace or commas; `comment` lists the characters that start a comment line (default `#`; `"#@"` for GROMACS xvg files) |
| `parse` | force a method: `regex`, `json`, `column`, `last_number` |

Fortran exponents (`1.5D-03`) are understood.

### `target` (required) and `relative`

The allowed error in the observable, absolute by default; `relative: true` makes it a
fraction of the observable's value at the starting setting.

### `cost` (optional)

How the cost of one run is measured. Default: wall-clock seconds. Otherwise a dict like
`observable` that reads a number from the output, for instance OpenFOAM's
`cost: {regex: "ExecutionTime\\s*=\\s*([0-9.]+)"}` or a solver that counts its own work.
Core-seconds (time × ranks) are the right cost for parallel runs; read the time and put
the rank count in the command if it varies.

### `failure` (optional)

When is a run not an answer? A non-zero exit code or a time-out always is a failure
(`allow_nonzero_exit: true` relaxes the first). In addition:

| key | |
|---|---|
| `failed_if_regex` | a pattern (or list) that marks a crash: `"FOAM FATAL"`, `"Traceback"` |
| `unconverged_if_regex` | a pattern that marks a run whose inner solver did not converge; the value is recorded but the run is never proposed or fitted unless `options.accept_unconverged` is true |
| `converged_if_regex` | a pattern that must appear for the run to count as converged: `"SIMPLE solution converged"`, `"convergence has been achieved"` |
| `path` | a file (or list of files) in the run directory to search as well |

The patterns are searched in stdout, stderr, the files the observable and the cost are
read from, and the files named in `path`. A solver whose output is redirected
(`simpleFoam > log.simpleFoam`) is therefore covered when the observable or the cost comes
from that log, or when `path: log.simpleFoam` is given.

### Files for each run

| key | |
|---|---|
| `template_dir` | a directory copied into every run directory (an OpenFOAM case, a LAMMPS deck). Every text file in it that contains `{{` is rendered |
| `files` | a mapping `source: destination` of single files to copy and render, the destination relative to the run directory |
| `render_files` | if given, only these files (relative to the run directory) are rendered |
| `workdir` | where run directories and the log go (default `errbudget_runs`, next to the job file) |
| `keep_runs` | keep run directories after the number is read (default true; set false to save disk) |
| `timeout` | seconds per run; a run that exceeds it is a failed run |
| `env` | environment variables for the command: `{OMP_NUM_THREADS: "1"}` |

### `options`

| key | |
|---|---|
| `configuration` | `v1.2.1` (default, the article's), `v1.2` or `v1.1` |
| `screen` | `true` runs the preliminary sweep: fewer runs, recommendations that keep frozen parameters at your value |
| `accept_unconverged` | `true` uses unconverged runs as answers (version 1.1 behaviour) |
| `floor_probe` | two or more settings (mappings of parameter values) refined beyond every candidate value; the spread of the observable over them is the round-off floor, and targets below ten times it are declined without searching. Recommended for single-precision codes |

Any other key is passed to the engine as an expert override (`resid_margin`, `joint_res`,
`max_iter`, …; see `docs/api.md`).

### `output`, `name`

`output` is the stem of the result files (`errbudget_result` → `errbudget_result.json`,
`errbudget_result.md`); `name` is a title for the report.

## The run log and resuming

Every run is appended to `<workdir>/runs.jsonl`: its parameters, value, cost, status
(`ok`, `unconverged`, `failed`), exit code, directory and wall time. Before running a
setting, the runner looks it up there, so

- `errbudget run` after an interruption continues where it stopped;
- `errbudget run` on a finished job answers everything from the log in seconds;
- changing the target or the candidate values and running again reuses every run that is
  still needed.

Every record carries a fingerprint of the command, the observable, cost and failure rules,
the environment and the template files. After you change any of those, the old records are
not reused (the log says how many it is ignoring); runs made with the same job file are.
`errbudget check` runs the starting setting again even when the log holds a failed attempt, so
a fixed command is actually tried. `errbudget run --fresh` ignores the log altogether; delete
`workdir` to start over completely. `errbudget runs job.yaml` prints the log as a table.

A time-out kills the whole process group of the command, so an `mpirun` that exceeded
`timeout` does not linger.

## Exit codes

`errbudget run` exits 0 when the target was met, 3 when it was not reachable with the values
given (the report says which parameter to extend), 1 on a job-file error.

## A complete example

```yaml
name: lid-driven cavity, OpenFOAM
command: blockMesh > log.blockMesh && simpleFoam > log.simpleFoam
template_dir: case                          # 0/, constant/, system/ with {{placeholders}}
parameters:
  - {name: n_cells,  kind: discretisation, values: [20, 40, 80, 160], start: 40}
  - {name: p_tol,    kind: nested-solver,  values: [1e-4, 1e-6, 1e-8], start: 1e-6}
  - {name: n_sweeps, kind: cost-only,      values: [1, 2, 3], start: 2}
observable:
  name: u_centre
  source: file
  path: postProcessing/probes/0/U          # lines like:  t (ux uy uz)
  column: 2                                # brackets are ignored when the line is split
target: 2e-3
relative: true
cost: {source: file, path: log.simpleFoam, regex: "ExecutionTime = ([0-9.]+) s"}
failure:
  failed_if_regex: "FOAM FATAL"
  converged_if_regex: "SIMPLE solution converged"   # searched in log.simpleFoam too (the cost's file)
env: {OMP_NUM_THREADS: "1"}
timeout: 7200
workdir: /scratch/me/cavity_runs
keep_runs: false
```
