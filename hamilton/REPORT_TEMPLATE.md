# Run report, cluster studies J1–J5

Copy this file, fill it in while the jobs run, and send it back with
`hamilton_results.tar.gz` (and `knobs_v6_results.tar.gz` once Section 6 of the runbook
has been done).

## Machine and versions

- Hamilton: node type (`sinfo`), OpenFOAM release and container tag (or module), Quantum
  ESPRESSO (`pw.x -v`), LAMMPS (`lmp -h | head -1`), GROMACS (`gmx_d --version`),
  Python / NumPy / SciPy (`python -c "import numpy, scipy, sys; print(sys.version, numpy.__version__, scipy.__version__)"`)
- Repository commit: `git rev-parse HEAD`; `git status --short` (should be clean apart from
  `config.sh`)
- `$WORK` and what it held before submitting (`ls $WORK`): anything moved aside, anything reused

## Smoke tests

```
(the five SMOKE TEST lines from */logs/smoke_*.out)
```

Anything changed to make them pass, as `git diff`.

## Jobs

- Final `bash status.sh` output.
- Core-hours per job (`sacct -X -u $USER -S <start date> --format=JobName,Elapsed,AllocCPUS,CPUTimeRAW,State`),
  wall time, number of resubmissions per job.
- Tasks that failed for good (two retries), each with the first lines of its `error.txt`.
- J2: Courant numbers of the spin-ups (`cat $WORK/dfg2d2_cases/mesh_*/SPUN_UP`); any
  `SPIN-UP FAILED` line; whether the 800 core-hour guard was reached and what remained.
- J1: `iters` and `converged` of `fo_na192`.
- J5: whether `$WORK/abc_hpc` held results before this run.

## Collect

- The `tail -n 3 results/j*.txt` lines.
- The checklist of Section 5 of the runbook: which items hold, which do not.

## Bringing the results back (if done on the article machine)

- `logs/v6/STATUS.txt` after `after` and `collect`.
- The line `frozen program verified ...` from B01.
- Whether the 2D-2h steps ran or were skipped.

## Deviations

Every change to code or configuration, as `git diff`; anything not run, and why.
