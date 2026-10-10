# Changelog

All notable changes to `knobs`. The article describes version 1.2; the comparisons
with earlier versions that used to be in the article are recorded here, with the
script that reproduces each.

## 1.2.0 — revision 3 of the article

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
