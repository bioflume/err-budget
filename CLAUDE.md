## Code style

- Write what an experienced research software engineer would keep after review: every line has a job.
- `src/errbudget/` is the record of the code that produced the results: do not edit it. Nor `PREREGISTRATION.md`, `paper/frozen_src/`, `hamilton/J5_abc3d/abc_hpc.py`, or the code of `hamilton/` and `studies/`, which ran as committed.
- Replay scripts live flat in `paper/` (the article cites them by path and line); shared logic in `paper/cases.py`, `pooled.py`, `baselines.py`, `paths.py`. Scripts read `data/`, write `data/*.json` and `figures/`, and run from `paper/` with `pip install -e .`. Import; never copy code between scripts or edit `sys.path`.
- A new configuration of the program is a line in `CONFIGS` of `paper/pooled.py`, not a new script. Every new result gets a step in `paper/replay_all.sh` and a row in `paper/README.md`.
- Comments explain why (equation, section of the article, numerical reason), never what. No banners, emoji or commented-out code.
- Docstrings: what the script computes, the section of the article it serves, how to run it, what it writes. No history; git keeps it.
- No `try/except` unless the handler does something specific. Never swallow an error or fall back silently.
- Printed output of a result script is a result: change it only on purpose, and regenerate `figures/<name>.txt`.
- Keep seeds and the released configuration `errbudget.V121`; unconverged runs count as failed.

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
