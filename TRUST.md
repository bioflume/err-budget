# Evidence for trusting errbudget

An AI coding agent wrote most of this repository (`AI_USE.md`). This file lists the evidence
that a reader can check without trusting the author or the agent, and what that evidence
found. It follows the rules of the UKRI project "Trustworthy Agentic Workflows in
Scientific Computing": the author owns the requirements and the acceptance tests, the checks
sit outside the agent's control, and failures are reported as findings.

Status: 26 November 2026, branch `paper`, engine 1.2.1 (`src/errbudget/_core.py`, SHA-256
55f0e994…).

## What a reader can run

| Check | Command | Expected |
|---|---|---|
| Acceptance tests | `pytest -q` | 52 passed |
| Specification of Section 4 | `python tests/test_spec_v121.py` | 14 passed |
| Human-owned files unchanged | `python tools/check_protected.py` | no difference from `PROTECTED.sha256` |
| Frozen program of the pre-registration | `python tools/check_frozen.py` | `paper/frozen_src` and the tag `v1.2.0-frozen` (commit 41d3db3) match the hash in `PREREGISTRATION.md` |
| Every number of the article | `bash paper/replay_all.sh`, then `git status paper/data paper/figures` | no stored result changes (`paper/REPRODUCTION_CHECK.md`) |

CI runs the first four on every push.

## The test designs, and which this code has

The project's work package 1 builds an acceptance rule from seven test designs. This table
says which of them this repository has.

| Design | Here |
|---|---|
| The community's tests | `tests/`: 52 tests; `test_claims.py` has one test per behavioural claim of the article, `test_spec_v121.py` one per rule of Section 4 |
| Conformity with established codes | the method runs five production codes it did not write (LAMMPS, GROMACS, ABINIT, Quantum ESPRESSO, OpenFOAM) against published or converged references (`THIRD_PARTY_RESULTS.md`) |
| Tolerances from conditioning and discretisation error | a target is graded only if it is at least ten times the uncertainty of its reference (`paper/cases.py`, `Case.admitted`) |
| Known answers, manufactured solutions | the development problems L2–L5 and the held-out problems H1 (manufactured solution) and H3 (exact solution) |
| New inputs from a generator | hidden checks, held outside this repository (below) |
| CPU–GPU agreement | not applicable: the engine runs on the CPU only |
| Per-result certificates | the engine reports an error estimate with each setting; it is an estimate with safety factors, not a rigorous bound, and it misses in the cases below |

## Pre-registered held-out tests

`PREREGISTRATION.md` fixed three held-out problems, the program, the targets and the grading
before any replay, and the article reports the results whatever they were (Table 5). The
frozen program still matches its registered hash (`tools/check_frozen.py`). The tag
`v1.2.0-frozen` had not been pushed to GitHub, so the CI check of it never ran before
26 November 2026.

The registrations of the three cluster studies (airfoil, reuse2d1, PRE comparison) are
commits of a project repository that is not available, so their dates cannot be checked
here. What remains is the dates written in the protocols, the dated amendments and notes of
the airfoil protocol (`studies/airfoil/notes/`), and the time stamps in the names of the
result packs (`DATA.md`), which all follow the registration dates. These are records the
author kept, not independent evidence.

## Hidden checks

A separate agent session, which did not write the engine, wrote a generator of synthetic
problems with known answers. It covers seven families: separable power-law errors,
pre-asymptotic coarse levels, a solver error that grows with mesh refinement, cancellation
between parameters, an inert parameter, round-off noise, and failed and unconverged runs.
The generator is kept outside this repository, so an agent working here cannot see it or
tune the engine to it. It grades with the rule of `PREREGISTRATION.md`.

Result for engine 1.2.1 (35 problems, seed 20261126, 97 test cases):

- **Searches.** 14,410 searches; 11,688 accepted, 2,610 declined although a setting met the
  target, 112 starts failed.
- **Misses.** 150 of the accepted searches missed. All 150 come from one test case
  (problem P08, pre-asymptotic, target 10⁻³), in which every start returns the same setting.
  Its true error is 1.23 times the target, while the engine estimated 0.85 times.
- **Elsewhere.** Every other family had no miss; the worst accepted error there was 0.82
  times the target.

## Planted faults

A separate agent (Claude Fable 5.2), which never saw the tests, planted 50 faults in the
engine and in `tune`: sign, normalisation and index errors, omitted terms, substituted
rules, and workflow faults. It also wrote 20 harmless variants, such as reordered sums. Each
was applied alone to a fresh copy of the code and run against the tests and the hidden
checks (21 problems).

| Check | Faults detected | Harmless variants flagged |
|---|---|---|
| Tests of this branch (`pytest`) | 36 of 50 (72%) | 0 of 20 |
| Hidden checks | 16 of 50 (32%) | 0 of 20 |
| Either | 38 of 50 (76%) | 0 of 20 |

The project's target is a design that detects at least 90% of faults with at most 5% false
alarms. The tests here do not meet it. Eleven faults pass both checks; the tests of `main`
miss the same ones, except F48, which does not apply there. They fall into three groups:

- **Inverted convergence flag.** Unconverged runs are used as answers by default (F48).
- **Weaker error estimates.**
  - The sign of the fitted coefficient is flipped (F05).
  - The residual of the fit is squared (F09).
  - The jointly refined change is left out of the lower bound (F28).
  - Under cancellation, the estimate takes the largest single term instead of the sum (F37).
- **Changes to the search or its cost accounting.** F08, F15, F23, F32, F41, F49.

The list, the harness and the results are kept with the hidden checks. On the branch `main`,
`tests/test_v122.py` (for the author to review, `AGENTS.md`) now catches six of these faults
(F05, F08, F09, F23, F37, F49); its tests were written knowing the faults, so this is not an
independent measurement.

## Independent review

Two reviewers, both Claude Sonnet 5.5 (a different model from the one that wrote the code),
read the code without the commit history, the change log or any test result. The agent
that prepared this file (Claude Opus 5.5) re-ran the findings marked "reproduced" below.

On the program that produced the results (`_core.py`, `api.py`) and the grading code:

- **The grading follows the pre-registered rule.** The reviewer re-implemented it from the
  text of `PREREGISTRATION.md` and got the stored H1 and H2 counts exactly.
- **An axis whose error cannot be measured is counted as error-free (reproduced).** On a
  ladder with three values whose middle run fails, a search started at the finest value
  certifies it. The true error there is 98 times the target. On 40 random tensors with
  three to five values per parameter and 12% of runs failed, 185 of 1,324 accepted
  settings missed, against none of 1,387 without failures; the reviewer did not finish
  counting how many of the 185 come from this cause. The hidden checks above use ladders of
  five or six values, so they do not reach this case.
- **The local-order check is skipped at the finest value of a ladder.** The pre-registered
  miss of Table 5 (H3, kinetic energy, 10⁻⁴, 216 of 216 searches) returns a setting at the
  finest value of both the mesh and the time step. With the check extended to the finest
  value, in a copy of the engine outside this repository, those 216 searches decline (no
  setting on the grid meets that target). The other H3 targets then accept fewer searches:
  96 of 216 instead of 216 for the kinetic energy at 10⁻³. The same change does not remove
  the hidden-check miss P08.
- **"No credit for cancellation" is narrower in the code than in the text.** Cancellation
  between parameters is credited unless the net error falls below 0.8 times the largest
  single term (`_core.py:431`); `PREREGISTRATION.md` says "no credit for cancellation".
- **`tune(..., relative=True)` scales the target by the value at the start setting,** which
  can be coarse, and runs the start setting twice. The article's replays do not use `tune`.

None of these changes a stored result. The engine on this branch stays as the article
describes it; the findings are fixed in engine 1.2.2 on the branch `main` (its default
configuration `V122`), whose configuration `V121` still runs this engine search for search.
The full reports are kept with the hidden checks.

## Code quality

Ruff (rules E, F, B, C90) and Radon on `src/errbudget/`:

- **Complexity.** Four functions of the engine have a cyclomatic complexity above 35:
  `refine_path` 44, `per_knob_error` 42, `joint_correct` 39, `run_protocol` 36; `tune` has 45.
- **Lint.** 94 warnings: long lines (46), lambda assignments (10), ambiguous names (9),
  complex functions (9), several statements on one line (8), unused imports (3), and nine
  "loop variable in a closure" (B023). The nine B023 warnings are callbacks used within the
  same iteration, so they do not change results.
- **Manual review.** Ten functions were drawn at random (seed 20261126) for the author's
  review.

## What this evidence does not show

- That the error estimate is a bound: it is not, and the hidden checks and Table 5 show
  misses.
- How the code was written: the agent's session logs are not in the repository
  (`AI_USE.md`).
- That the hidden checks are independent of the model that wrote the engine. The same model
  family (Claude) wrote both, in different sessions.
