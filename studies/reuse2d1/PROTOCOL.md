# Pre-registration: reuse and payback on a Reynolds-number study of 2D-1

Registered on 14 November 2026, before any run of this study. The airfoil study
(`../airfoil/`, amendment A6) ended without answering its questions on reuse and payback;
this study asks them on the article's OpenFOAM case 2D-1, where the program is known to
work at Re = 20. Changes after registration are amendments (§12), each committed alone
before the code or runs it affects.

## 1. Questions

- **R1 (reuse).** errbudget searches once, at the hardest member of a 13-member study in the
  Reynolds number, from the default start. Does the setting it returns meet the target at
  every member? Does it from every converged starting setting?
- **R2 (payback).** After how many production runs does the search pay for itself against
  the cautious setting, and does it within this study?
- **R3 (practice, information only).** How do the tutorial and the cautious settings fare at
  the same members, and a setting searched for at a member other than the hardest?

## 2. Problem

| | |
|---|---|
| Case | The article's 2D-1 case, unchanged (`paper/openfoam_cylinder.py` of the article repository, commit 4d32148): DFG cylinder in a channel, H = 0.41, D = 0.1, the twelve-block O-grid, parabolic inflow, `noSlip` walls and cylinder, p = 0 at the outlet; `simpleFoam` (SIMPLEC, `consistent yes`), U relaxation 0.9, `bounded Gauss linearUpwind` for U, `Gauss linear corrected` elsewhere; at most 2,000 iterations |
| Members | Re = Ū D / ν ∈ {10, 12.5, 15, …, 40} (13 members), Ū = 2 U_m / 3, ν = 10⁻³: the inflow maximum U_m = 0.015 Re is the only change from 2D-1 (U_m = 0.3, Re = 20). All are below the onset of vortex shedding |
| Code | OpenFOAM v2606 (ESI), Hamilton's module, driven through its input files only |
| Observable | the drag coefficient c_D = 2 F_x / (Ū² D · depth), from the written fields by the article's post-processing (pressure times face area plus the wall shear from the wall cells' centroids). c_L and the pressure difference are recorded, not graded |
| Cost | OpenFOAM's ExecutionTime of `simpleFoam` (CPU seconds), as in the article; cells × iterations is recorded too |
| Converged | `simpleFoam` reports "SIMPLE solution converged" before its iteration cap. A run that stops at the cap is unconverged; a run that fails is a failed run; neither is used as an answer |

## 3. Parameters

The article's lists and classes for 2D-1, unchanged:

| Parameter | Values | Class | h |
|---|---|---|---|
| `n_arc` (cells per quarter of the cylinder) | 8, 12, 16, 24, 32, 48, 64 | discretisation | 1/n_arc |
| `res_ctrl` (residualControl, p; U one tenth) | 10⁻², 10⁻³, …, 10⁻⁷ | nested solver | value |
| `lin_tol` (GAMG tolerance for p; U ten times it) | 10⁻⁴, …, 10⁻⁸ | nested solver | value |
| `n_coarsest` (GAMG nCellsInCoarsestLevel) | 10, 50, 250, 1000 | cost-only | 1/value |

The default start is the middle of each list: n_arc 24, res_ctrl 10⁻⁵, lin_tol 10⁻⁶,
n_coarsest 250. The 840 settings form the grid.

## 4. Runs

1. **The grid:** all 840 settings at every member, 10,920 runs.
2. **The references:** at every member, n_arc 96, 128 and 192 at res_ctrl 10⁻⁸, lin_tol
   10⁻¹⁰, n_coarsest 10, at most 20,000 iterations (39 runs).
3. **A check job** of three grid runs before the study (§11).

## 5. The reference and when a target is graded

At each member the reference value Q_ref is c_D at n_arc 192. From n_arc 96, 128 and 192 the
observed order p = ln(|Q128 − Q96| / |Q192 − Q128|) / ln 1.5 (ratio of the h steps 4/3 and
3/2 taken as 1.5 for both, a conservative simplification); if the three change
monotonically, p' = min(max(p, 1), 2), otherwise p' = 1. Its uncertainty is

u = |Q192 − Q128| / (1.5^p' − 1) / |Q192|.

A reference run that does not converge leaves its member without a reference. A target ε
(relative to Q_ref) is graded only if ε ≥ 10 u at every member that has a reference, and
only at the members that have one. At Re = 20, Q_ref is also compared with the benchmark
value 5.57953523384.

## 6. Targets

ε = 1 % (primary) and 0.3 % (secondary) of c_D.

## 7. The hardest member

At the default start, the relative change of c_D from n_arc 24 to n_arc 32 (the other
parameters at their default values) at each member; the member with the largest change is
the hardest. A member where either run did not converge is the hardest; among several, the
largest Re.

## 8. Procedures

- **(a) errbudget.** The released engine (version 1.2.1, `errbudget.tune`, the classes and h
  of §3), a relative target, from the default start, at the hardest member, replaying the
  stored runs (which gives the setting a live run would). Its setting is reused unchanged at
  every member. If it declines, the decline is the result.
- **(a, every start).** The same from every setting whose own run converged at the hardest
  member; each returned setting is reused at every member.
- **(b) tutorial.** n_arc 32, res_ctrl 10⁻², lin_tol 10⁻⁶, n_coarsest 10 (the article's
  tutorial tolerances at n_arc 32).
- **(c) cautious.** The finest value of every list: n_arc 64, res_ctrl 10⁻⁷, lin_tol 10⁻⁸,
  n_coarsest 10.
- **(a′) another member (information).** (a) searched from the default start at each other
  member, reused at every member.

A setting meets the target at a member if its run there converged and
|Q − Q_ref| ≤ ε |Q_ref|.

## 9. Outcomes

- **R1.** For (a): the members where it meets each graded target and the largest ratio of
  error to ε. For (a, every start): the share of (start, member) pairs that meet the target,
  the declines and the distinct settings. For (b), (c) and (a′): the same, as information.
- **R2.** The search cost S is the CPU cost of the distinct runs the search of (a) made
  (failed and unconverged runs at their own cost). K = S / cost(a) at the hardest member,
  µ = cost(c) / cost(a) there, and n* = K / (µ − 1), as in the article. For the study:
  N* = S / mean over members of (cost(c) − cost(a)), the number of members after which the
  search has paid for itself; it pays within this study if N* ≤ 13. If (a) is not cheaper
  than (c) at most members, the search does not pay against (c), and that is the result.
- CPU time varies between Hamilton nodes; the same quantities are also given in cells ×
  iterations.

## 10. Outcomes reported as negative

- (a) misses a graded target at any member, with its cause.
- errbudget declines from the default start at the hardest member.
- A target not graded (§5), with the member whose reference carries u.
- N* > 13 against (c) at a target: the search does not pay within this study.

## 11. Check job, budget and stop rules

The check job runs three grid runs (Re 20 at n_arc 8 and 12, Re 40 at n_arc 8, res_ctrl
10⁻⁴, lin_tol 10⁻⁶, n_coarsest 10) and passes if all three converge and c_D at Re 20,
n_arc 12, is within 3 % of the benchmark value. Predicted cost: 40–90 core-hours
(`README.md`); guard 200 core-hours, past which the study stops and what remains is
reported.

## 12. Amendments

None.
