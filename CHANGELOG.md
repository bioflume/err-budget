# Changelog

## 2.1.0 — engine 1.2.2 (27 November 2026)

Engine 1.2.2 corrects three rules of 1.2.1 that an independent review found (`TRUST.md`),
and is the new default. The configuration `"v1.2.1"` (`errbudget.V121`) still runs the
engine exactly as the article describes it: on 86,001 searches over synthetic problems, every
proposal, decision, run count, cost and error estimate is bit for bit that of 2.0.0. The
article's results are on the branch `paper` and do not change.

- **Engine 1.2.2: an axis whose error cannot be measured is not error-free.** An axis with
  fewer than three successful runs and no finer run (for example a three-value ladder whose
  middle run failed, searched from its finest value) was left out of the error estimate, and
  a setting with 98 times the target error was certified. It is now charged the change to
  its nearest finite coarser run (more where the values are closely spaced: that change
  extrapolated to first order), and the search declines (`unmeasured_axis`) if the axis has
  no finite change at all (`strict_unfit`).
- **Engine 1.2.2: the local-order check runs at the finest value of a ladder** (from the three
  values ending there), where 1.2.1 skipped it (`local_at_finest`). On the pre-registered
  held-out problem H3 this turns the 216 searches that 1.2.1 accepted and that missed the
  kinetic-energy target 10⁻⁴ (no setting on the grid meets it) into declines.
- **Engine 1.2.2: no credit for any cancellation between parameters.** With
  `cancel_rule="none"`, 1.2.1 summed the errors in magnitude only when their signed sum fell
  below 0.8 times the largest term; errors of +1, +1 and -1 were credited two thirds of
  their size. 1.2.2 sums them in magnitude whenever their signs differ (`cancel_strict`).
- **`errbudget.V122`** is V121 with these three switches on; `tune` and job files use it
  unless told otherwise (`configuration="v1.2.1"`, `options.configuration: v1.2.1`).
- **What it changes on the article's stored runs** (every graded test case and the held-out
  problems H1-H3, every start): V121 accepts 28,257 of 36,315 searches and misses 216 (all
  on H3); V122 accepts 23,971 and misses none. The price is caution: on the development
  problems V122 accepts 22,625 searches against 25,869, and on H1-H3 1,346 against 2,388.
  Most of the difference comes from the local-order check at the finest value, the
  correction that removes the H3 misses (`TRUST.md`).
- **`tune(configuration=...)`** names the configuration: "v1.2.2" (default), "v1.2.1", "v1.2",
  "v1.1". Keyword overrides are added to it. Before, `configuration: v1.2` in a job file and
  `tune(..., **V12)` ran with the joint-check resolution of 1.2.1 (1, not 3), because `V12` has
  no `joint_res` key; `full_configuration(name)` now lists every option of a configuration.
- **Fixed: `tune(relative=True)`** scaled the target by the value at the start setting, so a
  coarse start certified a miss of the relative target. The search is now repeated with the
  target scaled by the returned value less its estimated error, until it holds. The start
  setting was also run twice; every setting now runs once. A failed or zero start value is an
  error with a message instead of a target of `nan`.
- **Fixed: the floor probe.** A failed probe run is left out instead of stopping the tuning
  with a traceback; fewer than two usable probe runs, or a probe setting that does not give
  every parameter, is an error with a message.
- **Fixed: a declined result says why** (`TuneResult.declined`, in the report and in Markdown):
  a failed start, a floor, an unmeasured or unbounded parameter, or no setting within the
  values given.
- **Fixed: job files are checked more strictly.** Refused, with a message: unknown keys in a
  parameter, `observable` or `cost` (a misspelt `knd:` was read as a discretisation
  parameter, `regx:` as the last number of the output); `relative` and `keep_runs` that are not
  true or false (`"false"` was read as true); a non-finite target; `column` below 1 (0 read the
  last column); `parse` without the key it needs; a list for `converged_if_regex` (it crashed
  every run); `render_files` that is not a list; `options` or `env` that is not a mapping;
  floor-probe settings that do not give every parameter; parameters named `run_dir`, `run_id`
  or `job_dir`; a field given twice; invalid YAML, JSON or TOML (a traceback before).
- **Fixed: runs of a command.** An observable that reads as `nan` or `inf` is a failed run
  (it was logged as ok); a `render_files` entry missing from the run directory fails the run
  (the file was left unrendered); a run that timed out is run again when the job's `timeout`
  is longer; a failed run is charged its recorded cost (`RunFailed(cost=...)`); where `files`
  copies each file is part of the fingerprint of the run log.
- **Fixed: reports** print parameter values in full in Markdown (they were rounded to six
  digits) and list frozen parameters in the terminal report too.
- **Not changed:** the search ranks its moves as before; the order of `os.walk` in the
  fingerprint and numbers such as `"1e3"` read as floats are as in 2.0.0.

- **Removed from job files.** The field `solver` (`module:function`, `file.py:function`); from
  Python, call `tune`. A job file that still has it is refused as an unknown field. Run logs
  written by 2.0.0 are still reused.
- **Removed kind aliases.** `d`, `n`, `c` and `solver`; a job file that uses one is refused with
  a message that lists the accepted spellings. The alias `discretization` is now documented.
- **Removed `python -m errbudget.cli`.** Run the command as `errbudget` or `python -m errbudget`.
- **Test problems.** `AdvectionDiffusion`, which nothing used, is removed, and `make_synthetic`
  moved from `errbudget.problems` into `tests/test_synthetic.py`, its only user. Also gone: the
  problems' `runner()` methods and the private module `_asymptotic`, which only used each
  other; `AdvectionDiffusion(precond=...)`, which had no effect; the extra `examples`.
- **Fixed: the floor probe of `tune`** reads its runs like every other run. Without
  `observable=` it took the first value of a returned dict, which could be the `"converged"`
  flag, and then declined nothing. A probe run that did not converge is now left out.
- **Fixed: `failure:`** given as a string or a list is refused with the job-file message
  instead of a Python error; the exit code is still 1.
- **Fixed: `target: true`** in a job file is refused; it was read as an absolute target of 1.0.
- **Tests.** Redundant, always-passing and never-running tests are gone; the suite takes under
  a minute.
- **Planned, after the article is accepted.** Remove the engine's dead code
  (about 75 lines: `refine_path(n_levels >= 2)`, `cancel_rule="path"`,
  `resid_mode="additive"`, `unmeasured_rule`, `shares()`, `model_abs_error()`,
  `propose(mode=...)`, `run_protocol(mode=...)`) and the ablation-only options `walk_mode`,
  `phase2`, `screen_at` and `audit_factor` (about 65 lines). The results will be the same, and
  "the engine is unchanged" will then read "the engine gives the same results".

## 2.0.0 — the package is now called errbudget

Only names change. The engine is the same: `_core.py` is byte for byte that of 1.3.0, and
`errbudget.ENGINE_VERSION` is still 1.2.1, the version the article describes.

- **Package and command.** `pip install errbudget`, `import errbudget`, and the command
  `errbudget` (`init`, `check`, `run`, `report`, `runs`, `version`; also `python -m errbudget`).
- **`Knob` is now `Parameter`**, the second argument of `tune` is called `params`, and
  `JobSpec.knobs()` is now `JobSpec.to_parameters()`. Error messages say "parameter".
- **Defaults that carried the old name.** `errbudget init` writes `errbudget.yaml`; runs go to
  `errbudget_runs` and the results to `errbudget_result.json` and `.md`; a saved result
  records `errbudget_version`. A job file that sets `workdir` and `output` is unaffected. To
  resume a tuning that 1.3.0 started with the default run directory, set
  `workdir: knobs_runs`.
- **Repositories.** <https://github.com/bioflume/errbudget>; the article's repository is
  <https://github.com/bioflume/errbudget-paper>.
- Kept: the engine's own names (for example `LocalBudget.per_knob_error`) and the job-file
  format, which already used `parameters`.
- **Project pages.** A logo; a documentation site built from `docs/` and published at
  <https://bioflume.github.io/errbudget/>; a contributing guide, a code of conduct, a security
  policy, issue forms and a pull-request template. CI also runs on Python 3.14, and the
  `dev` extra installs SymPy, which some of the test problems need.

The entries below keep the old name.

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
repository: <https://github.com/bioflume/errbudget-paper/blob/main/CHANGELOG.md>.
