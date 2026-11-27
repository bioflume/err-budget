# Rules for coding agents

Coding agents (Claude Code and others) wrote most of this repository under the direction of
its author, Gökberk Kabacaoğlu, who is responsible for every change and every result
(`AI_USE.md`). These rules apply to any agent working here. They follow the rules of the
UKRI project "Trustworthy Agentic Workflows in Scientific Computing": humans own the
requirements and the acceptance tests, agents write code, and the checks stay outside the
agents' control.

## Files you must not change

The files matched by `PROTECTED` are the requirements, the acceptance tests and the engine
that produced the article's results. Do not edit, move, delete or regenerate them, and do
not edit `PROTECTED`, `PROTECTED.sha256`, `tools/check_protected.py`, `AGENTS.md`,
`AI_USE.md`, `TRUST.md`, `.github/CODEOWNERS` or `.github/workflows/`. CI fails when a protected file
differs from `PROTECTED.sha256`, and only the author updates that manifest.

- The engine `src/errbudget/_core.py` is engine 1.2.2. Its configuration `V121` must keep
  reproducing the article's engine 1.2.1 run for run: every correction is an option that
  `V121` leaves off (CONTRIBUTING.md, "The engine rule"). Its hash is in `PROTECTED.sha256`,
  so a change to it needs the author's approval.
- `tests/test_spec_v121.py` and `tests/test_claims.py` pin the rules of the article, and
  `tests/test_v122.py` the corrections of engine 1.2.2. If one fails, the code is wrong or
  the rule is wrong; never change the test to make it pass.

## When a check fails or a task cannot be done honestly

Stop and report. Say which check failed, with its output, and what you think the cause is.
Do not weaken a test, special-case an input, catch and hide an error, loosen a tolerance or
skip a check to make a task pass. A task you cannot complete within these rules is a
finding, not a failure: report it.

New tests you write go in a new file that the author reviews and, if accepted, adds to
`PROTECTED`. A test an agent wrote does not count as acceptance evidence until then.

## What to record

- Every commit you make ends with a `Co-Authored-By:` line that names the model and its
  version.
- When you run anything that produces a number for the article or the documentation, say
  which script, commit and environment produced it.
- Report results as they are: failures and negative results included.

## Code style

- Write what an experienced research software engineer would keep after review: every line
  has a job.
- The package is `src/errbudget/`: `api.py` (`tune`, `Parameter`), `jobfile.py` and
  `runner.py` (job files and runs of a command), `cli.py` (`errbudget init/check/run`),
  `report.py`, `problems/` (the test problems). Examples in `examples/` import the package;
  they never edit `sys.path` or copy code.
- Validate input where it enters (`Parameter`, job files, the command line) with a clear
  message; not on internal calls.
- No `try/except` unless the handler does something specific (report a failed run, retry a
  parse). Never swallow an error or fall back silently.
- Comments explain why (equation, rule of the article, numerical reason), never what. No
  banners, emoji or commented-out code.
- Numpydoc docstrings on the public functions (`tune`, `Parameter`, `run_protocol`,
  `load_job`, `run_job`), with what each returns; one line or none on private helpers.
- No logging configuration or prints inside the package, except the command-line tool and
  `verbose=True`.
- A new option of the program is a key of the configuration (`V122`) or of the job file,
  documented in `docs/`, with a test in `tests/`.
- Run `pytest tests -q` and `python tools/check_protected.py` before every commit.
