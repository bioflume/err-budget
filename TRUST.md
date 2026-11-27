# Evidence for trusting errbudget

An AI coding agent wrote most of this repository (`AI_USE.md`). This file lists the evidence
that a reader can check without trusting the author or the agent, and what that evidence
found. It follows the rules of the UKRI project "Trustworthy Agentic Workflows in
Scientific Computing": the author owns the requirements and the acceptance tests, the checks
sit outside the agent's control, and failures are reported as findings.

Status: 27 November 2026, branch `main`, errbudget 2.1.0 (release), engine 1.2.2. The article
describes engine 1.2.1; the configuration `v1.2.1` (`errbudget.V121`) runs it unchanged
(below). The branch `paper` holds the article's own record and its `TRUST.md`.

## What a reader can run

| Check | Command | Expected |
|---|---|---|
| Acceptance and interface tests | `pytest -q` | 114 passed |
| Human-owned files unchanged | `python tools/check_protected.py` | no difference from `PROTECTED.sha256` |
| The command-line example | `cd examples/04_command_line && errbudget check job.yaml && errbudget run job.yaml --quiet` | runs to the end |

CI runs all three on every push. `tests/test_v122.py` (39 of the 114 tests) was written by
the agent for engine 1.2.2; the author reviewed and accepted it on 27 November 2026, and it
is in `PROTECTED`.

## The article's engine is unchanged under `v1.2.1`

Each correction of engine 1.2.2 is a switch that `V121` leaves off. Two comparisons check
that `V121` still runs the article's engine:

- **Synthetic problems.** The engine 1.2.1 file and the 1.2.2 file were loaded side by side
  and run on the same problems from every start: 86,001 searches over 653 test cases, with
  short ladders and 12% failed runs among them. None differed: not the setting,
  the decision, the decline reason, the number of runs, the cost or the error estimate.
- **The article's stored runs.** Every graded test case and the held-out problems H1–H3 were
  replayed with engine 1.2.2 under `V121`. Of the 83 rows, 80 have a stored counterpart in the
  article's tables (`paper/data/pooled_v1.2.1.json`, `heldout_replay_current.json` on the
  branch `paper`), and all 80 equal it.

## What engine 1.2.2 changes

On the article's stored runs (36,315 searches: every graded test case and H1–H3, every
start), graded as the article grades them:

| Configuration | Development problems: accepted / missed | Held-out H1–H3: accepted / missed |
|---|---|---|
| v1.2.1 | 25,869 / 0 of 33,333 | 2,388 / 216 of 2,982 |
| v1.2.2 | 22,625 / 0 | 1,346 / 0 |
| v1.2.2 without the local-order check at the finest value | 25,074 / 0 | 2,324 / 216 |
| v1.2.2 without the strict short-axis rule | 23,343 / 0 | 1,406 / 0 |
| v1.2.2 without strict cancellation | 22,813 / 0 | 1,358 / 0 |

- **Misses.** Engine 1.2.2 removes the 216 misses of 1.2.1, all of them on H3 (kinetic
  energy, target 10⁻⁴, which no setting on the grid meets): those searches now decline.
- **Cost of the caution.** 1.2.2 accepts 13% fewer searches on the development problems and
  44% fewer on H1–H3. A declined search is safe but leaves the user without a setting.
- **Which correction causes it.** Most of the difference comes from the local-order check at
  the finest value, the correction that removes the H3 misses.

A first version of the short-axis rule scaled a ten-fold step down by 9. It added two misses
on H2 that 1.2.1 did not have, and was corrected before release: the rule never charges less
than 1.2.1.

## Hidden checks of the engine

A separate agent session, which did not write the engine, wrote a generator of synthetic
problems with known answers. The generator is kept outside this repository and grades with
the article's pre-registered rule. On 35 problems (97 test cases, 14,410 searches):

| Configuration | Accepted | Missed | Worst accepted error |
|---|---|---|---|
| v1.2.1 | 11,688 | 150 | 1.23 × target |
| v1.2.2 | 10,926 | 150 | 1.23 × target |

All 150 misses, under both configurations, come from one pre-asymptotic test case in which
every start returns the same setting. There the mesh error changes sign between the coarse
and the fine values, and the observed order overstates the true one. Engine 1.2.2 does not
fix this case; no correction was designed for it.

## Planted faults

A separate agent (Claude Fable 5.2), which never saw the tests, planted faults in the code
and wrote harmless variants. Each was applied alone to a fresh copy and run against the
tests and the hidden checks (21 problems). A check detects a fault when it fails; it raises
a false alarm when it fails on a harmless variant.

**Round 1, engine 1.2.1 (50 faults, 20 variants; 49 and 19 apply to this branch).** The
tests of 2.0.0 detected 37 of 49 (76%), with no false alarm. The new tests of
`tests/test_v122.py` were written afterwards, knowing these faults. They catch six of the
ten that escaped (F05, F08, F09, F23, F37, F49). Still escaping:

- F15, F32 and F41, which change how the search ranks its moves;
- F28, which leaves the jointly refined change out of the estimate (a term that rarely
  matters).

**Round 2, engine 1.2.2 (50 new faults, 20 new variants, written blind to the tests,
including 15 in the job-file and runner layer).**

| Check | Faults detected | Harmless variants flagged |
|---|---|---|
| Tests of this branch (`pytest`, 114 tests) | 25 of 50 (50%) | 0 of 20 |
| Hidden checks (engine, 21 problems, V122) | 2 of 50 (4%) | 0 of 20 |
| Either | 25 of 50 (50%) | 0 of 20 |

By file, the tests caught 15 of 27 faults in the engine, 6 of 9 in `runner.py`, 2 of 8 in
`api.py` and 2 of 6 in `jobfile.py`. Several of the 25 escapes are in code added for 1.2.2:

- the rescaling of relative targets (G07, G22, G40, G43);
- the short-axis rule (G03, G11, G19);
- duplicate keys in JSON job files (G30).

Others, such as G47, are older paths. In G47, floats rounded in the run-log key make the
tolerances 10⁻⁸ and 10⁻¹⁰ share a run. The tests of this branch are the weakest evidence
here.

The project's target is at least 90% of faults detected with at most 5% false alarms.

## Independent review

Two reviewers, both Claude Sonnet 5.5 (a different model from the one that wrote the code),
read the code without the commit history, the change log or any test result. The agent
that prepared this file (Claude Opus 5.5) re-ran their main findings. All confirmed
findings are fixed in 2.1.0 (`CHANGELOG.md`):

- **Engine.**
  - An axis whose error cannot be measured counted as error-free (98 times the target
    error, certified).
  - The local-order check was skipped at the finest value of a ladder.
  - "No credit for cancellation" credited cancellation unless the net error fell below 0.8
    times the largest term.
- **`tune`.**
  - `configuration: v1.2` ran with the joint-check resolution of 1.2.1.
  - A relative target was scaled by the value at the start, and the start setting ran twice.
  - A failed start gave a target of `nan`.
  - A failed floor-probe run ended the search with a traceback.
  - A declined result did not say why.
- **Job files and runs.**
  - Misspelt keys were ignored.
  - `relative: "false"` read as true.
  - `column: 0` read the last column.
  - A list for `converged_if_regex` crashed every run.
  - Unchecked `render_files` left files unrendered.
  - A `nan` observable was logged as a good run.
  - A timed-out run stayed cached after the time limit was raised.
  - Failed runs lost their cost.
  - Duplicate keys and invalid YAML were not reported as job-file errors.

Not changed, from the reviewers' unconfirmed suspicions:

- the order of `os.walk` in the run-log fingerprint (a missed reuse, not a wrong one);
- numbers such as `"1e3"` being read as floats;
- run directories left after an interrupt;
- the quoting of a command given as a list.

The full reports are kept with the hidden checks.

## Code quality

Ruff (rules E, F, B, C90) and Radon on `src/errbudget/`, 2.0.0 → 2.1.0:

- **Complexity.** Cyclomatic complexity of the largest functions:
  - `JobSpec.from_dict` 112 → 142;
  - `tune` 39 → 64;
  - `per_knob_error` 42 → 60;
  - `refine_path` 44 → 45;
  - `joint_correct` 39 → 41;
  - `run_protocol` 36 → 36.

  The corrections made the validation and `tune` larger. They are candidates for splitting
  into smaller functions.
- **Lint.** 183 → 213 warnings, all of style classes (162 long lines). The eight "loop
  variable in a closure" warnings (B023) are callbacks used within the same iteration.

## What this evidence does not show

- That the error estimate is a bound: it is not. The hidden checks still show 150 misses in
  one pre-asymptotic case under 1.2.2.
- That 1.2.2 is the better default for every user. It removes the article's held-out misses
  at the price of more declines, and the author decides the default.
- How the code was written: the agent's session logs were not kept (`AI_USE.md`).
- That the hidden checks and the planted faults are independent of the model family that
  wrote the engine: Claude models wrote all of them, in separate sessions.
