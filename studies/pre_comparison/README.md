# errbudget against probabilistic Richardson extrapolation (PRE) on OpenFOAM 2D-1

A head-to-head comparison on the stored runs of the article's 2D-1 study. `PROTOCOL.md` fixes the design
and was committed before any PRE result was graded. `REPORT.md` gives the results.

| Step | Command | Output |
|---|---|---|
| Fetch the authors' PRE code and its environment | `sh fetch_pre.sh` | `ext/`, `.venv/` |
| V0: the authors' unit tests | `cd ext/sparse-…; PYTHONPATH=src ../../.venv/bin/python -m unittest test.test_sparse_pre` | `validation/unittest.log` |
| V1: the published cubature example (rerun in `validation/rerun/`) | `cd validation/rerun; PYTHONPATH=../../ext/sparse-…/src ../../.venv/bin/python ../../ext/sparse-…/src/run_cubature_experiments.py all`, then `validate_spre.py published` | `validation/published.log`, `results/validation.json` |
| V2: `spre_white.py` against the authors' code | `PYTHONPATH=ext/sparse-…/src .venv/bin/python validate_spre.py agree` | `validation/agree.log`, `results/validation.json` |
| Arm 1 | `.venv/bin/python arm1.py` | `results/arm1_rows.csv` |
| errbudget against Table S-3.1, errbudget-3p | `.venv/bin/python check_replay.py` | `results/errbudget_check.json` |
| Arm 2 | `.venv/bin/python arm2.py` | `results/arm2_rows.csv` |
| Exploratory kernels | `PYTHONPATH=ext/sparse-…/src .venv/bin/python explore_kernels.py` | `results/explore_kernels.json` |
| Tables and criteria | `.venv/bin/python summarize.py` | `results/summary.json`, `results/tables.md` |
| Figure | `.venv/bin/python make_figure.py` | `fig_pre_2d1.pdf` |

The study and the program are read from the article's repository (`ERRBUDGET_REPO`, default
`~/OpenResearch/knobs-paper/knobs-paper`, commit 4d32148). Nothing there is modified.
