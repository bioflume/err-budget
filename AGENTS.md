# Rules for coding agents

Coding agents (Claude Code and others) wrote most of this repository under the direction of
its author, Gökberk Kabacaoğlu, who is responsible for every change and every result
(`AI_USE.md`). These rules apply to any agent working here. They follow the rules of the
UKRI project "Trustworthy Agentic Workflows in Scientific Computing": humans own the
requirements and the acceptance tests, agents write code, and the checks stay outside the
agents' control.

## Files you must not change

The files matched by `PROTECTED` are the requirements, the acceptance tests, the program
that produced the results and the record of the results: `src/errbudget/`, the
pre-registration and the frozen program, the stored runs and printed outputs
(`paper/data/`, `paper/figures/`, `hamilton/results/`, `studies/*/summary/`), the study
protocols and their notes. Do not edit, move, delete or regenerate them, and do not edit
`PROTECTED`, `PROTECTED.sha256`, `tools/`, `AGENTS.md`, `AI_USE.md`, `TRUST.md`,
`.github/CODEOWNERS` or `.github/workflows/`. CI fails when a protected file differs from
`PROTECTED.sha256`, and only the author updates that manifest.

- `src/errbudget/` is the record of the code that produced the results. The engine
  `_core.py` is byte-identical to the one on `main`.
- `tests/test_spec_v121.py` and `tests/test_claims.py` pin the rules of the article. If one
  fails, the code is wrong or the article is wrong; never change the test to make it pass.
- A replay writes its results over the stored ones on purpose, so that `git diff` shows any
  difference. If a replay changes a stored result, report the difference; do not commit it.

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

- Write what an experienced research software engineer would keep after review: every line has a job.
- The code of `hamilton/` and `studies/` ran on the cluster as committed: do not edit it.
- Replay scripts live flat in `paper/` (the article cites them by path and line); shared logic in `paper/cases.py`, `pooled.py`, `baselines.py`, `paths.py`. Scripts read `data/`, write `data/*.json` and `figures/`, and run from `paper/` with `pip install -e .`. Import; never copy code between scripts or edit `sys.path`.
- A new configuration of the program is a line in `CONFIGS` of `paper/pooled.py`, not a new script. Every new result gets a step in `paper/replay_all.sh` and a row in `paper/README.md`.
- Comments explain why (equation, section of the article, numerical reason), never what. No banners, emoji or commented-out code.
- Docstrings: what the script computes, the section of the article it serves, how to run it, what it writes. No history; git keeps it.
- No `try/except` unless the handler does something specific. Never swallow an error or fall back silently.
- Printed output of a result script is a result: change it only on purpose, and regenerate `figures/<name>.txt`.
- Keep seeds and the released configuration `errbudget.V121`; unconverged runs count as failed.
- Run `pytest -q`, `python tools/check_protected.py` and `python tools/check_frozen.py`
  before every commit.

Job scripts (Hamilton, Durham) start from this and add only what the job needs:

    #!/bin/bash
    #SBATCH --partition=shared
    #SBATCH --time=<hh:mm:ss>
    #SBATCH --cpus-per-task=<n>
    module purge
    module load openfoam/v2606
    source $HOME/venvs/airfoil/bin/activate
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    python driver.py "$@"

No echo banners, timing, `hostname`, `cd $SLURM_SUBMIT_DIR` or input checks: `sacct` and `seff` record job metadata.
