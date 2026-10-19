# Changelog

## 1.3.0 — tooling around the article's engine

The engine (`knobs/_core.py`, the defaults `V121`) is byte-for-byte version 1.2.1, the
version the article describes; `knobs.ENGINE_VERSION` records it. New around it:

- **The `knobs` command**: `init` (annotated job file), `check` (validate and run the
  starting setting once), `run`, `report`, `runs`, `version`.
- **Job files** (YAML, JSON or TOML) for solvers that run as a shell command: placeholders
  in the command and in template files, per-run directories, the observable and the cost
  read from stdout or a file (regex, JSON key, table column, last number), failure and
  convergence patterns, time-outs, a run log that resumes an interrupted tuning. Python
  functions can be named in a job file too (`solver: module:function`).
- `TuneResult.to_dict / to_json / load / markdown`.
- Documentation: user guide, job-file reference, recipes for OpenFOAM, LAMMPS, GROMACS,
  Quantum ESPRESSO/ABINIT and MPI codes, API reference; a command-line example.
- NumPy 2 allowed (`numpy>=1.20`); Python 3.9+.

## 1.2.1 — the article's version (revision 3)

- The joint-refinement check includes every parameter whose one-step change exceeds the
  resolution (`joint_res=1.0`) instead of three times it. Found by the pre-registered
  held-out problem H2; registered as an amendment before the third held-out problem ran.
  `knobs.V12` reproduces the frozen version 1.2.

## 1.2.0 — frozen for the held-out tests

- Joint-refinement check, local-order check, no credit for cancellation, unconverged runs
  treated as failed, resolution floor, residual-gated safety factor, full mode by default.

## 1.1 — revision 2 of the article

- `unfit="order1"`, the `noise_floor` check, `models=[..."exponential"...]`.

The complete history with the comparisons behind each change is in the article's
repository: <https://github.com/bioflume/knobs-paper/blob/main/CHANGELOG.md>.
