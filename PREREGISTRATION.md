# Pre-registration: held-out tests of knobs 1.2

**Purpose.** Referee point 3 of revision 3: "freeze version 1.1, choose one or two
problems or codes that played no part in development, and report the results
whatever they are."

**Changes since the referee's request.** Version 1.2 differs from 1.1 in the changes
listed in `CHANGELOG.md`, all of which were designed on the development problems
(L2–L5, OpenFOAM 2D-1 and 2D-2, LAMMPS, GROMACS, ABINIT and Quantum ESPRESSO). Those
changes are:

- the joint-refinement check;
- the local-order check;
- no credit for cancellation;
- unconverged runs treated as failed;
- the resolution floor;
- the residual-gated safety factor as the default;
- the full mode as the default.

For that reason the frozen program is 1.2, not 1.1.

**What this document fixes.** This document was committed to the repository before
any held-out run was replayed. It fixes:

- the program;
- the configuration;
- the problems, their parameters and declarations;
- the targets;
- the reference values;
- the grading rule.

Any later change to any of these would be visible in the repository history.

## The frozen program

- **Code:** `src/knobs/` at the git tag `v1.2.0-frozen`. The SHA-256 of the concatenated SHA-256 sums of `src/knobs/*.py` is given at the end of this file, and `paper/exp62_heldout.py` checks it before replaying.
- **Configuration:** the default of version 1.2, `knobs.V12`:
  - residual-gated safety factor (α = 5);
  - joint-refinement check;
  - local-order check;
  - resolution of one tenth of the target;
  - no credit for cancellation between parameters.
- **Failed runs:** a run whose nested solver did not meet its own convergence criterion is treated as failed.
- **Mode:** full mode (no preliminary screening).

## Problems

None of these problems, codes or observables was used while the program was developed.

### H1 — 3D variable-coefficient Poisson problem, scikit-fem (local)

**Module:** `paper/heldout_fem.py`.

**Problem.**
- −∇·((1+x)∇u) = f in the unit cube.
- Manufactured solution u = sin πx sin πy sin πz + xyz.
- P2 tetrahedra.
- Solver: conjugate gradients preconditioned by smoothed-aggregation AMG (pyamg).

**Observables** (exact values):

| Observable | Exact value |
|---|---|
| mean of u | 8/π³ + 1/8 |
| flux of (1+x) ∂u/∂x through the face x = 1 | 2(1/4 − 4/π) |

**Parameters** (refinement scalar *h* in parentheses):

| Parameter | Values | Declared class | Model |
|---|---|---|---|
| n | 3, 4, 6, 8, 12, 16, 24 (*h* = 1/n) | discretisation | power |
| cg_rtol | 1e-2 … 1e-10 (*h* = value) | nested-solver | power |
| amg_theta | 0, 0.08, 0.25 | cost-only | — |

**Cost:** deterministic work count (CG iterations × nonzeros, plus the AMG hierarchy).

**Note on the solver at n = 24 (before the freeze, before any replay).** The runs at
n = 24 exceeded the memory of the development machine. From that size on, the forms
are assembled over chunks of elements, and the L2 projection of the boundary data is
solved by conjugate gradients to a relative residual of 1e-14 instead of directly.
At n = 8 this changes both observables by less than 1e-14 relative
(`python paper/heldout_fem.py check 8`).

**Targets** (relative):

| Observable | Targets |
|---|---|
| mean | 1e-3, 3e-4 |
| flux | 3e-2, 1e-2, 3e-3 |

### H2 — lid-driven cavity at Re = 1000, OpenFOAM simpleFoam v1912 (local)

**Module:** `paper/heldout_cavity.py`.

**Observables:** extrema of the centreline velocities. Reference values are from
the benchmark of Botella & Peyret (1998), given to seven digits:

| Observable | Reference value |
|---|---|
| u_min | −0.3885698 |
| v_max | 0.3769447 |
| v_min | −0.5270771 |

**Parameters:**

| Parameter | Values | Declared class |
|---|---|---|
| N | 24, 32, 48, 64, 96, 128, 192 (*h* = 1/N) | discretisation |
| res_ctrl | 1e-3 … 1e-7 | nested-solver |
| lin_tol | 1e-4, 1e-6, 1e-8 | nested-solver |

**Convergence:** a run that does not meet `res_ctrl` within 30,000 iterations is
reported as not converged.

**Cost:** cells × outer iterations.

**Targets** (relative): 3e-2, 1e-2, 3e-3 for each observable.

### H3 — decaying 3D ABC flow, OpenFOAM pimpleFoam, current release (Hamilton)

**Module:** `hamilton/J5_abc3d/abc_hpc.py`.

**Problem.** u = (sin z + cos y, sin x + cos z, sin y + cos x) e^(−νt) in (0, 2π)³,
with all boundaries cyclic, ν = 0.1, T = 2.

**Observables** (exact values at T):

| Observable | Exact value |
|---|---|
| ek | 1.5 e^(−2νT) |
| eps | 3ν e^(−2νT) |
| umax | √6 e^(−νT) |

**Parameters:**

| Parameter | Values | Declared class |
|---|---|---|
| N | 16, 24, 32, 48, 64, 96 | discretisation |
| dt | 0.04, 0.02, 0.01, 0.005 | discretisation |
| p_tol | 1e-3, 1e-5, 1e-7 | nested-solver |
| n_outer | 1, 2, 3 (*h* = 1/n) | nested-solver |

**Cost:** core-seconds.

**Targets** (relative):

| Observable | Targets |
|---|---|
| ek | 1e-3, 1e-4 |
| eps | 3e-2, 1e-2, 3e-3 |
| umax | 1e-3, 3e-4 |

## Grading (identical to the pooled table of the article, `paper/pooled.py`)

**Replay.** For every observable and target, the program is replayed from every
usable start, i.e. every setting whose run produced a value.

**Counts.** Each replay is counted as one of:

| Count | Meaning |
|---|---|
| start failed | the user's own run failed or did not converge |
| accepted | the program returned a setting it accepts |
| missed | accepted, but the true relative error of the returned setting exceeds the target |
| declined, reachable | declined, although some converged setting on the grid meets the target |
| declined, unreachable | declined, and no converged setting on the grid meets the target |

**Penalty.** For each accepted setting, the penalty is its cost relative to the
cheapest converged setting that meets the target. The median and the worst are
reported.

**Reporting.** The results are reported in full in the article, whatever they are,
together with version 1.1 on the same tensors for comparison. No problem, target or
observable listed here will be dropped.

## Frozen code hash

(filled in at the time of the freeze, see the git tag)
