# Changelog

All notable changes to `errbudget`, which was called `knobs` until revision 5 of the
article. The article describes version 1.2.1; the comparisons with earlier versions that
used to be in the article are recorded here, with the script that reproduces each. The
entries below the rename keep the old name, as it was at the time.

## Results of the cluster reruns (November 2026; program unchanged: 1.2.1)

- **Refined 2D-2 study (J2).** `hamilton/J2_dfg2d2/analyse_j2.py` grades against the run at
  $n_\mathrm{arc}=128$, $\Delta t=1.25\times10^{-4}$: the planned reference at 192 diverged
  (each of its 32 time steps reached the limit of 50 outer correctors). Its uncertainty, from
  the mesh change 96 to 128, admits eight cells. The case `2D-2h` (1,620 settings from 180
  runs; 1,260 with a value, 1,233 converged) replaces the coarser study in Tables 3, 4 and 7
  (`exp70_tables.py`, `ERRBUDGET_TABLE3_2D2=refined`, the default); `data/pooled_*.json` hold
  both studies. In the eight admitted cells v1.2.1 accepts 8,259 of 9,864 searches and misses
  none. Below the admission threshold it accepts, from all 72 starts for the drag at
  $3\times10^{-3}$ and all 64 for the lift at $10^{-2}$, settings at $n_\mathrm{arc}=32$ whose
  finer run diverged, and every one misses (by 1.44 and up to 2.39 times the target): the
  charge for a failed finer run assumes that the later changes halve (condition (A1) of
  Proposition 1). The engine is not changed; treating such a parameter as unmeasured unless
  three values through the current one ran removes these misses in replay and is left to a
  later version.
- **Totals with the refined study:** 56 cells, 32,091 replays, 24,842 accepted starts (305
  distinct settings), no accepted setting outside an admitted target.
- **J4.** The observables do not depend on the number of ranks (LAMMPS 1-64 MPI ranks,
  GROMACS double precision 1-16 thread-MPI ranks); the timings grew with the number of ranks
  and are not used.

## Cluster jobs: J2 inlet fix and reruns (November 2026; program unchanged: 1.2.1)

Only the cluster drivers changed (`hamilton/`, and failure messages in `paper/qe_silicon.py`,
`paper/gromacs_water.py`, `paper/exp55_parallel.py`); `src/errbudget` is unchanged.

- **Container launches and local paths.** `config.sh` binds the node-local `/local` (where
  SLURM puts `$TMPDIR`, and Open MPI its session directory) into the OpenFOAM container, and
  reads `config.local.sh` (paths of codes in their own environments; the old `KNOBS_*`
  names are accepted).
- **Slowed-down runs.** mpirun runs OpenFOAM with `--bind-to none`, and J2 stops a run
  whose linear solves take more than 0.1 s each over 15 minutes (`ofmodern.run`,
  `ofmodern.TooSlow`: retried and resumed, not a failed setting). A resumed J2 run is
  charged from the checkpoint it resumed from, so redone work is not counted twice; a
  segment that fails before its first time step is retried, not stored as a failed setting.
- **J4.** Jobs load `gcc/11.2` (the C++ runtime of the GROMACS build); whole-node tasks ask
  for 64 tasks; each LAMMPS task checks one cheap setting first and reports LAMMPS's
  messages (`exp55_parallel.LAST_ERROR`); `lmp_run` rejects a run whose log reports fewer
  processes than ranks; `submit.sh` refuses results without a value.
- **Resubmission.** `farm left` counts unfinished tasks; `submit.sh` skips meshes (J2) and
  studies (J4) that are complete; the J2 smoke test runs on one node, exercises the watchdog
  and fails if a run failed.
- **Collection.** `JOBS="2 4" bash collect_job.sh` runs only the analyses of the reruns;
  results moved aside (`before_*`, `failed_*`) are left out of `raw.tar.gz`.

- **J2: analytic inlet on mapped meshes.** `mapFields -consistent` also maps fixedValue
  patches, so every mesh initialised from a coarser one ran with the coarser mesh's inlet
  profile (flow rate off by up to 0.4 %). `spin()` restores and checks the analytic inlet
  after mapping and records it; meshes and results from before the fix are refused.
- **J2: spin-ups.** One SLURM job per mesh with a one-day limit; runs per mesh after their
  spin-up; a launch failure (no time step) no longer halves the spin-up time step;
  `SPUN_UP` records the mesh actually mapped from; `analyse_j2.py` reports whether the
  reference runs converged and refuses an unconverged main reference.
- **Task farm and monitoring.** Stale `error.txt` kept as `error_prev.txt`, attempts logged,
  started-but-unfinished tasks counted; `coreh.sh` and `status.sh` report core-hours.
- **J3, J4.** A run with no output raises (retried, with the program's error message); the
  smoke tests check finite results.

## Rename: `knobs` is now `errbudget` (program unchanged: 1.2.1)

The program's name describes what it returns, an error budget. Only names changed:

- **Package and imports.** `src/knobs` moved to `src/errbudget`, and every script, test,
  example and cluster job imports `errbudget` (`errbudget.V121`, `errbudget.tune`, ...).
  `_core.py`, the engine, is byte-identical; elsewhere in the package only docstrings, the
  header line of the printed report (`errbudget: recommendation for ...`) and the names
  below changed.
- **`Knob` is now `Parameter`**, the second argument of `tune` is called `params`, and the
  error messages say "parameter". The examples, the tests and `paper/exp23_mlhybrid.py` use
  the new name; the frozen program keeps `Knob`.
- **Environment variables** `KNOBS_*` are now `ERRBUDGET_*` (for example `ERRBUDGET_JOBS`,
  `ERRBUDGET_ROOT`, `ERRBUDGET_SMOKE`). On Hamilton, `config.sh` now defaults to the
  micromamba environment `errbudget` and the scratch directory `/nobackup/$USER/errbudget`;
  set `ERRBUDGET_ENV=knobs` or `WORK=/nobackup/$USER/knobs` to keep using an environment or
  directory made before the rename.
- **Repositories.** This repository is `errbudget-paper`, the library `errbudget`
  (<https://github.com/bioflume/errbudget-paper>, <https://github.com/bioflume/errbudget>).
- **Kept under the old name**, because they record the program as registered:
  `PREREGISTRATION.md`, the tags `v1.2.0-frozen` and `v1.2.1`, `paper/frozen_src/knobs`,
  and the registered J5 driver `hamilton/J5_abc3d/abc_hpc.py`, which still reads
  `KNOBS_SMOKE` (the J5 smoke test sets both names). `paper/exp62_heldout.py` loads the
  frozen program under its registered name, checks its hash against `PREREGISTRATION.md` as
  before, and registers it as `errbudget` for the replay scripts.
- **Metadata.** `.zenodo.json` now gives version 1.2.1 and the affiliation of the article
  (Department of Computer Science).
- **Checked:** all 52 tests pass; `exp70_tables.py`
  reproduces every table of revision 5 byte for byte; `exp62_heldout.py` verifies the frozen
  program and writes the same `data/heldout_replay*.json`, with the frozen program and with
  `--current`. With `--current` it prints a new source hash, because names and docstrings in
  the package changed; this hash was never registered, and the frozen one is unchanged.

## Revision 5 of the article (program unchanged: 1.2.1)

The program did not change. The article and its scripts did:

- **A proposition on when the estimate bounds the error** (Section 4.8 of the article).
  Every charge of the budget is the sum of a geometric series that continues a measured or
  fitted change (Lemma 1). The estimate bounds the error when, along every parameter, the
  changes beyond the current value decay at least as fast as one of these series, the
  changes of all parameters have one common sign, and the non-additive part of the error is
  at most (F_s - 1) times the sum of the errors along the parameters (Proposition 1). With
  exact separable power laws, the program accepts if and only if the finest setting meets
  0.8 eps and otherwise declines as `not_reached` (Corollary 1).
  `paper/exp81_bound_check.py` (step A21c) checks the three statements against the code and
  runs the synthetic counterexamples of the Remark: same-sign higher-order terms (missed by
  up to 8%) and a bilinear interaction with a parameter at its finest value (up to 28%).
- **The article reports the released program only.** The comparisons with earlier versions
  moved to Section S15 of the supplement and to the sections below: the version-1.1 column
  of Table 4, the version rows of Table 5 (the frozen program stays, as the joint-check
  threshold 3 rho), the version-1.1 column of Table 6, the version-1.1 columns of Table A.1,
  and the version comparisons in the text. `exp70_tables.py` writes the supplement tables
  `tab_versions_summary.tex`, `tab_versions.tex` and `tab_estimate_v11.tex`.
- **Results restated in terms of the components of the program.** `exp71_dfg2_detail.py`
  adds v1.2.1 with cancellation credited (2D-2 drag at 3e-3: 116 of 172 accepted starts
  miss, all at one setting, where a mesh error of +1.46 eps and a time-step error of
  -1.20 eps are credited as cancelling while the true error is +1.03 eps) and prints the
  budget there. `exp80_v6_numbers.py` adds the ABINIT and Quantum ESPRESSO budgets of the
  fits alone, without the local-order check, and (part 5) the budget at the two L5 settings
  that miss with a constant safety factor: at the finest mesh, where the local-order check
  has nothing to check, the fits underestimate the error by a factor of 1.35, and only the
  gated factor lifts the estimate (0.997 eps with the constant factor, 1.11 eps with the
  gated one, true error 1.08 eps).
- **Section 9 and the program summary** describe the package (<https://github.com/bioflume/knobs>,
  1.3.0, engine 1.2.1, the `knobs` command and job files) and this repository with its three
  reproduction tiers. The 2D-2 grid takes about 40 CPU-hours (revision 4 said 24).
- `paper/README.md`: steps A21b and A21c, and the map of every result of revision 5 to its
  script. `reproduce.sh`: the running time of the 2D-2 grid.

## Revision 4 of the article (program unchanged: 1.2.1)

The program did not change. The repository did:

- **New analyses** (`paper/run_v6_analyses.sh`, steps A01–A23, B01–B07, C01) for the
  referee report on version 5: the ABINIT reference uncertainty with the Quantum ESPRESSO
  recipe (`exp75`), the PPPM order mechanism (`exp74`), the accounting of Table 2 against
  Table 3 (`exp76`), relTol in 2D-1 (`exp77`), the GROMACS floor with v1.2.1 (`exp78`),
  unconverged runs start by start (`exp79`), the collection of every result
  (`exp82`), ablations from v1.2.1 as well as from v1.2 (`pooled.py`, `exp70_tables.py`),
  and `distinct` (accepted settings) in every pooled and held-out row.
- **Cluster jobs** (`hamilton/`): the J1 function object now writes (it was configured to
  write every 20,000th field write and never did); J2 spin-ups at Courant 0.65–0.87 with
  retries at half the time step; refusal to reuse results of the old drivers; `collect_job.sh`.
  Details in `hamilton/README.md`.
- **Results of Part V** (everything that needs no cluster) regenerated on 2026-10-19 from
  the stored tensors with NumPy 2.5.3 and SciPy 1.18.1: every count matches the revision-3
  files. `abinit_reference_v6.json` and `pppm_mechanism_v6.json` need ABINIT and LAMMPS and
  are not yet in `paper/data`; the cluster results (J1–J5, held-out problem H3) are not yet
  in `hamilton/results`.
- **Layout.** `hamilton/README.md` is now the complete runbook for the cluster studies, with
  `hamilton/REPORT_TEMPLATE.md`; `paper/README.md` describes the three reproduction tiers;
  the library for users of the method moved to its own repository
  (<https://github.com/bioflume/knobs>).

## 1.2.1 — revision 3 of the article (released default)

- **The joint-refinement check now includes every parameter whose one-step change exceeds
  the resolution** (`joint_res=1.0`), instead of three times the resolution.
  - Found by the pre-registered held-out problem H2 (OpenFOAM cavity at Re = 1000). There,
    for v_max at 3e-2, version 1.2 accepted a setting from 9 starts that missed by up to 3%,
    because both one-step changes were below 0.3 ε and the check was skipped.
  - Registered as an amendment in `PREREGISTRATION.md` before the third held-out problem was
    run.
  - Default configuration: `knobs.V121`. `knobs.V12` reproduces the frozen version
    (git tag `v1.2.0-frozen`, copy in `paper/frozen_src`).
- Effect on the development problems (`paper/pooled.py v1.2.1`; supplement, Section S15):
  - 17,610 starts accepted instead of 17,612; the only cell whose number of acceptances
    changes is the GROMACS mixed-precision energy at 1e-4 (69 instead of 71);
  - no misses in either version, and the median of the cell medians of the penalty stays
    at 2.26;
  - 12 of 53 cells change their returned settings. The largest changes are on L2, where
    the larger joint check measures interactions that stop the search at finer settings:
    the median penalty of the dissipation at 3e-3 rises from 1.22 to 2.16, and the worst
    penalty of the energy at 7e-3 from 15.6 to 61.9.
- Effect on the held-out problem H2 (not a held-out result, because the amendment was made
  after the replay): the nine misses of v_max at 3e-2 disappear without a change of
  acceptances; the median penalty of u_min at 3e-2 rises from 1.00 to 1.69.

## 1.2.0 — revision 3 of the article (frozen for the held-out tests)

The default configuration is `knobs.V12`. Every number in the article is for it unless
stated otherwise. `knobs.V11` reproduces version 1.1 exactly, so that the comparisons
below can be re-run.

### Added

- **Joint-refinement check** (`joint_check=True`).
  - Before a setting is accepted, every parameter that is not cost-only and has a resolved one-step change is refined together by one position. This costs one run.
  - The program compares the change with the sum of the one-at-a-time changes. The difference *I* is the interaction term, and it is added to the estimate.
  - The check was proposed by the referee of revision 2. It removes the misses that came from interacting errors of opposite sign (DFG 2D-2, synthetic test `tests/test_v12.py::test_interacting_errors_of_opposite_sign`).
- **Local-order check** (`local_order_check=True`).
  - Applies along every discretisation parameter. It uses three consecutive values through the current one, or the current value and the two finer ones at the coarse end, or the three values ending at the current one when the finer run failed.
  - It gives the observed order and extrapolates the remaining error with it. Changes that do not decrease mean the parameter's error is unbounded, so the program refines it or declines.
  - This catches candidate lists that start before the asymptotic range (ABINIT at 8 Ha) and observables whose error stalls between grids (L5).
- **No credit for cancellation** between parameters (`cancel_rule="none"`). If the signed sum of the per-parameter errors is smaller than 0.8 times the largest term, the magnitudes are summed instead.
- **Unconverged runs are failed runs.** A solver may return `{"...": value, "converged": False}` or raise `knobs.RunFailed`. Such a run is recorded with its cost, never proposed and never fitted. `tune(..., accept_unconverged=True)` restores the old behaviour.
- **Resolution floor** (`res_frac=0.1`). Changes below a tenth of the target, or below the round-off floor if one was measured, are treated as zero by the checks.
- **Declines are reported more precisely.** `declined` is one of `start_failed`, `unbounded_axis`, `unmeasured_axis`, `not_reached` or `below_floor`.

### Changed

- **The residual-gated safety factor is now the default**: `resid_margin=5.0`, i.e. F_s = 1.25 (1 + 5 q̄).
- **Full mode is the default** (`screen=False`). The preliminary sweep is still available as `screen=True`.
  - In revision 2 the sweep was the default.
  - It explains the misdeclaration result the referee queried. With the sweep, a parameter declared as nested solver is never frozen, and so it is refined to the cheapest valid setting. Without the sweep, every declaration gives the same recommendation (supplement, Table S2).

### Measured effect (pooled over all development cases, `paper/pooled.py`)

See Tables S9 and S10 of the supplement (revision 5) and `paper/data/pooled_*.json`.
Version 1.1, replayed in the full mode on the same stored runs:

- accepted 19,888 starts and missed from 456 of them, which returned seven distinct
  settings in five cells: the y-coordinate of the L5 minimum at 3e-4 (348), the ABINIT
  energy at 1e-2 and 3e-3 (48 and 36, all from the starts at 8 Ha), the Quantum ESPRESSO
  energy at 1e-4 (23) and LAMMPS at 8,000 ions and 1e-6 (1); 21 of its accepted settings
  were unconverged runs; median penalty 1.54;
- the L5 and ABINIT misses are the two symptoms the local-order check detects (an error
  that stalls between grids, a list that starts before the asymptotic range); the Quantum
  ESPRESSO misses come from an exponential fit that extrapolates too far, and the gated
  safety factor alone removes them (v1.1+gated); the LAMMPS miss is caught by more than one
  of the new checks;
- versions 1.2 and 1.2.1 remove every miss at the cost of 2,278 acceptances and a median
  penalty about 1.5 times higher (2.26).

Elsewhere (supplement, Section S15):

- held-out H1 and H2 with version 1.1: every start of H1 accepted at four targets without
  a miss; on H2, with unconverged runs accepted, 105 of 105 starts accepted in eight cells,
  with misses for v_max at 3e-2 (18 starts, one setting) and u_min at 3e-2 (3);
- OpenFOAM 2D-2 below the admission threshold: for the drag at 3e-3, version 1.1 missed
  from 168 of 195 accepted starts (median penalty 0.5) and version 1.2 with cancellation
  credited from 115; for the lift at 1e-2, version 1.1 missed from 196 of 198;
- GROMACS mixed precision without the floor probe: version 1.1 missed from 6 of 9 accepted
  starts at 1e-5 (up to 1.58 eps) and from 6 of 6 at 3e-6 (up to 3.2 eps);
- ABINIT: the budget of version 1.1 at 12 Ha and 3^3 k-points gives 99.8% of the energy
  error to the cutoff and 99.7% of the force error to the k-points; with the power-law
  model, version 1.1 declined every energy target and accepted the force only from the 48
  starts at 24 and 28 Ha with the two coarsest k-meshes;
- the estimate of version 1.1 was within the tightest graded target at 29 settings whose
  true error exceeded it, all for the y-coordinate of the L5 minimum (Table S11).

## 1.1.0 — revision 2

- **A parameter that cannot be fitted** (fewer than three successful runs along it) is charged its one-step change, extrapolated at first order (Eq. 4 of the article).
  - In version 1.0 such a parameter entered only through the one-step bound.
  - It removed:
    - the optimistic estimates at the coarse end of the candidate lists (spacing study: safety at a spacing of 1.2 rose from 90% to 100%);
    - the OpenFOAM 2D-1 misses of version 1.0 (17, 90 and 56 starts in three cases, up to 2.7 ε).
- **A parameter whose finer run failed** is charged its last measured change rather than treated as converged. This removed the two LAMMPS misses at 262,144 ions and 10⁻⁶, which were caused by the PPPM grid exceeding the memory limit.
- **Round-off floor probe** (`floor_probe=`): targets below ten times the measured spread are declined. This removed the misses of the mixed-precision GROMACS build (6 of 9 and 6 of 6 accepted settings at 10⁻⁵ and 3 × 10⁻⁶).
- **Exponential error model** for plane-wave cutoffs (`Knob(..., model="exponential")`). With the power law, version 1.0 declined every ABINIT energy target.
- **Residual-gated safety factor** (`resid_margin`), optional in this version.

## 1.0.0 — revision 1

- Offset power-law fit per parameter.
- Signed error budget.
- One-step bound.
- Constant safety factor 1.25.
- Monotone two-phase search.
- Preliminary sweep with frozen-parameter check.
- Cost-only sweep at the answer.

## 0.x

Development versions; see `DESIGN_NOTES.md`.
