# Changelog

All notable changes to `knobs`. The article describes version 1.2.1; the comparisons
with earlier versions that used to be in the article are recorded here, with the
script that reproduces each.

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
- Effect on the development problems (`paper/pooled.py v1.2.1`):
  - 17,610 starts accepted instead of 17,612;
  - no misses in either version;
  - 12 of 53 cells change.

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

See Table 3 of the article and `paper/data/pooled_*.json`.

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
