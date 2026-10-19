# Cluster studies on Hamilton: the runbook

Everything in this folder runs on Durham's Hamilton cluster (SLURM). Nothing else in the
repository needs the cluster. The five jobs here produce the results of the article that
need more than the 2-core machine used for the rest of it, and this file is the complete
procedure: setting up, checking, submitting, watching, collecting, and bringing the
results back into the article's analyses.

Read it top to bottom once. The commands are meant to be pasted in order.

| Job | What it answers | Core-hours | Wall time |
|---|---|---|---|
| `J1_dfg2d1` | 2D-1: forces from OpenFOAM's own function object against the article's post-processing; drag convergence to n_arc 192; the wall distance of the two cell centres; the whole 840-setting grid with a current OpenFOAM release, with the function object and the stopping reason of every linear solve; relTol 0 and 0.01 at n_arc 32 | ~20 | < 1 day |
| `J2_dfg2d2` | 2D-2 refined to n_arc 96 (reference at 128 and 192), with the averaging window and the relaxation time as parameters, forces from the function object | 300–600 | 1–3 days |
| `J3_qe` | the Quantum ESPRESSO grid with a current release (no `snprintf` shim) | ~3 | hours |
| `J4_parallel` | LAMMPS at 262,144 ions and GROMACS at 13,824 waters on 1–64 ranks: does communication change the cost ordering or the recommendation? | ~15 node-hours | < 1 day |
| `J5_abc3d` | **held-out 3D test**: decaying ABC flow with an exact solution, run with the frozen program as pre-registered in `../PREREGISTRATION.md` | ~60 | 1 day |

The estimates come from per-cell timings on the article's machine (about 1.4 µs per cell
and outer iteration for `pimpleFoam`); Hamilton's cores are faster.

## What is in this folder

```
config.sh            paths, partitions, account, OpenFOAM mode -- the only file to edit
setup_env.sh         one-off: micromamba environment, OpenFOAM container, GROMACS (double)
run_all.sh           submits J1-J5 (or their smoke tests with "smoke")
status.sh            tasks done / with errors, per job
collect.sh           runs the five analyses and packs hamilton_results.tar.gz
collect_job.sh       runs collect.sh as a batch job (the J1 and J2 analyses take 1-2 hours)
common/              farm.py (restartable SLURM task farm), ofmodern.py (OpenFOAM driver
                     for current releases), submit_array.sh
J1_dfg2d1/ ...       one folder per job: the driver (*_hpc.py), submit.sh, analyse_j*.py
J3_qe/pseudo/        the pseudopotential for J3
REPORT_TEMPLATE.md   what to write down while the jobs run and send back with the results
```

Each job sits in its own folder, runs as SLURM arrays of independent, restartable tasks,
writes one small `result.json` per task under `$WORK`, and has an analysis script that
turns those results into the numbers and tables of the article.

## 0. Rules that apply throughout

1. **The held-out test (J5) stays as registered.** Do not change `src/knobs/`,
   `../PREREGISTRATION.md`, `J5_abc3d/abc_hpc.py` or the H3 targets. J5 runs once. Its
   results are reported as they come out, good or bad.
2. **Do not change numerical settings**: grids, candidate lists, targets, tolerances, time
   steps, meshes, references or grading rules. What may be changed is what stops a job from
   running: paths, module or container names, partitions, account, memory and time limits.
   Record every change (`git diff`).
3. **A failed run is a result.** Never edit, delete or regenerate a result file to make a
   failure go away. A diverged or unconverged setting is data. An `error.txt` beside a task
   means an infrastructure problem: fix the cause and resubmit.
4. **One commit, traceable.** Work on a commit of this repository (`git rev-parse HEAD`
   goes into the report), so that every result can be traced to the code that produced it.
5. **Stop and ask instead of improvising** when a smoke test fails for a reason other than
   configuration, when the frozen-program check fails, when J2 passes 800 core-hours, or when
   the same task fails twice for the same reason.

## 1. One-off setup (login node, 30–60 minutes)

```bash
git clone <this repository> knobs-paper && cd knobs-paper/hamilton
nano config.sh        # WORK (scratch, default /nobackup/$USER/knobs); SB_EXTRA (account, if
                      # yours needs one); OF_MODE (container or module, see below)
bash setup_env.sh     # safe to re-run: finished steps are skipped
```

`setup_env.sh` installs micromamba in `~/bin`; creates the conda environment `knobs`
(Python 3.11, NumPy, SciPy, SymPy, pytest, a current Quantum ESPRESSO, LAMMPS with OpenMPI,
scikit-fem, pyamg, and this package in editable mode); pulls the OpenCFD container
`opencfd/openfoam-default:2406` into `$WORK`; and builds GROMACS 2024 in double precision
(the conda-forge GROMACS builds are single precision, whose round-off floor is too high
for the targets of J4).

**OpenFOAM: container or module.** The default is the container. To use Hamilton's own
module instead, set `OF_MODE=module` and `OF_MODULE=` to the name that
`module avail openfoam` shows; the drivers then run in the module's environment.

**Earlier runs under the same `$WORK`.** The task farm skips every task that already has a
`result.json`, so results of an earlier run would be reused. Check `ls $WORK` before
submitting:

- `dfg2d1_hpc`, `dfg2d2_hpc`, `dfg2d2_cases` from before revision 4: move them aside
  (`mv $WORK/dfg2d1_hpc $WORK/old_rev3/`, and so on); do not delete them. Their drivers
  changed, and the J1 and J2 submit scripts refuse to run over results of the old drivers
  and print the `mv` command.
- `qe_hpc`, `parallel_hpc`: keep them; their drivers did not change and finished tasks are
  reused.
- `abc_hpc` (J5): if it already holds results, **J5 has already run**. Do not delete it and
  do not run it again: report it and collect what is there.
- `smoke`: may stay; the J1 and J2 smoke tests clear their own smoke folders first.

## 2. Smoke tests (test partition, 15 minutes each; all five must pass before Section 3)

```bash
bash run_all.sh smoke                     # or: (cd J1_dfg2d1 && bash submit.sh smoke) per job
# when they have run:
grep -h "SMOKE TEST" */logs/smoke_*.out   # five lines, all PASSED
```

What each one checks:

- **J1** runs two function-object levels, two grid runs and one relTol run on tiny meshes.
  It passes only if OpenFOAM's `forceCoeffs` output is read as a number in `fo_na8` and in
  `grid_0_1_2_0`, a grid run has stopping-criterion counts, and the relTol run wrote
  `residuals.npz`. The log also prints the J1 analysis tables for the smoke tasks.
- **J2** checks the spin-up, `mapFields`, `decomposePar` and a parallel run in the container.
  It does not run the analysis (the smoke windows are too short for a shedding period).
- **J3, J4, J5** run one or two tiny tasks through the whole pipeline.

If J1 fails, read `J1_dfg2d1/logs/smoke_*.out` and any `$WORK/smoke/dfg2d1_hpc/*/error.txt`.
The likely causes are the container, or the `forceCoeffs` keywords of your OpenFOAM
release: the function object is written by `common/ofmodern.py` (`force_coeffs_fo`), and
its output is read by `read_coefficients`, which looks the columns up by name (`Cd`, `Cl`).
A keyword fix there is allowed; record it.

| Smoke test fails because of | Do |
|---|---|
| a path, module, partition, account or container name | fix `config.sh`, re-run, record the change |
| the driver or the function-object reader | the smallest fix that lets it run, without changing numerical settings; record the diff; if none exists, stop and report |

## 3. Submit

```bash
bash run_all.sh           # J1-J5; prints the job ids
```

J2 submits one 32-core spin-up job first (coarse to fine meshes); its run arrays depend on
it and start when it ends, even if a mesh level failed. Every task writes `result.json` when
it finishes.

## 4. Monitor and resubmit (daily, for 1–3 days)

```bash
squeue -u $USER
bash status.sh            # tasks done / with errors, per job; lists error.txt files
sacct -X -u $USER -S <start date> --format=JobName,Elapsed,AllocCPUS,CPUTimeRAW,State
```

- **Time-outs.** Run `bash submit.sh` again in that job's folder, as often as needed.
  Finished tasks are skipped and the others resume from their checkpoints.
- **`error.txt`.** Read it, fix the infrastructure cause, resubmit. After two failed retries
  of the same task, stop retrying it and report it.
- **J2 spin-up.** Check `J2_dfg2d2/logs/spin_*.out`. A failed level is retried at half the
  time step, twice; after a time-out, resubmitting resumes at the time step it was using
  (kept in `SPIN_ATTEMPT`). A level that still fails prints `SPIN-UP FAILED for n_arc N`
  and its runs report "not spun up": do not change time steps by hand, report it.
  `$WORK/dfg2d2_cases/mesh_*/SPUN_UP` records the largest Courant number of each spin-up;
  expect 0.4–0.9 and report any value above 1.
- **J2 cost guard.** Runs at the coarsest time steps on the finest meshes have Courant
  numbers of about 2–10 and may hit the limit of 50 outer iterations in many steps. That is
  expected and part of the study. Check the accumulated core-hours daily with `sacct`; if
  J2 passes **800 core-hours**, stop submitting J2 and report which tasks remain.
- **J1 `fo_na192`** may stop at the 20,000-iteration limit without meeting residualControl
  1e-8. Report `iters` and `converged`; do not raise the limit.
- A setting that diverges or does not converge needs nothing: it is recorded as a result.

## 5. Collect (when `status.sh` shows every job complete, or only reported failures left)

```bash
bash collect_job.sh       # batch job; the J1 and J2 analyses take 1-2 hours on one core
# when it ends:
ls -la hamilton_results.tar.gz && tail -n 3 results/j*.txt
```

`collect.sh` writes `results/j1.txt` … `results/j5.txt`, the tensors and JSON files the
article's analyses read (`results/*.pkl`, `results/*.json`), the raw per-task files
(`results/raw.tar.gz`, under 150 MB), and packs them into `hamilton_results.tar.gz`.

Check before sending it back:

- `results/j1.txt`: numbers (not `nan`) in the `Cd (FO)` column for all ten levels; the
  wall-distance table; the stopping-criterion table; the relTol table; replay lines for
  three tensor variants ("v1912 (article)", "current, fields", "current, function object").
- `results/j2.txt`: the accounting line `settings (base run x relaxation x window): 1620`;
  no `NOT FINISHED` list; the reference and its uncertainty; the Courant table; the line
  `admitted targets ...`; replay lines for five configurations, with non-admitted targets
  marked `[not admitted]`.
- `results/j3.txt`: 480/480 settings and replay lines for v1.2.1 and v1.1.
- `results/j4.txt`, `results/j5.txt`: no traceback.

## 6. Bring the results into the article (on the machine with the article's analyses)

Copy `hamilton_results.tar.gz` to the machine that holds this repository with the stored
tensors (`paper/data/`), keeping the original on Hamilton. Then, from `paper/`:

```bash
cd paper
bash run_v6_analyses.sh check                               # Python stack, solvers, git tag
bash run_v6_analyses.sh after /path/to/hamilton_results.tar.gz   # 2-3 hours
bash run_v6_analyses.sh collect                             # RESULTS_SUMMARY.md + archive
```

`after` runs `exp60_hamilton.py` (unpacks into `hamilton/results`, prints J1–J5, replays
H1–H3 with the frozen program and with v1.2.1, regenerates the tables; the line
`frozen program verified` must appear), then merges the refined 2D-2 study (`2D-2h`) into
every pooled file, the baselines, the accounting, the unconverged counts, the decline
reasons and the final tables (`exp70_tables.py`; Table 3 uses the refined study by default,
`KNOBS_TABLE3_2D2=refined`). If `hamilton/results/j2_tensor.pkl` lacks the case
specification, the 2D-2h steps are skipped and reported: re-run the patched
`J2_dfg2d2/analyse_j2.py` on Hamilton. `collect` writes `RESULTS_SUMMARY.md` in the
repository root and `knobs_v6_results.tar.gz` beside it. The steps that do not need the
cluster (`run_v6_analyses.sh before`) are described in `../paper/README.md`; run them first.

## 7. What to send back

1. `hamilton_results.tar.gz` (from Section 5) and, after Section 6, `knobs_v6_results.tar.gz`.
2. A report following `REPORT_TEMPLATE.md`: machine and versions, the five smoke-test lines,
   the final `status.sh` output, core-hours per job, resubmissions, tasks that failed for
   good with the first lines of their `error.txt`, the Courant numbers of the J2 spin-ups,
   and every change to code or configuration as `git diff`.

## Notes

- **Held-out test (J5).** J5 is the held-out test of the frozen program (version 1.2) and,
  by the amendment of 11 October 2026, of version 1.2.1. Its problem, parameters, targets
  and grading rule are fixed in `../PREREGISTRATION.md`, together with the SHA-256 of the
  frozen source.
- **Costs.** Costs are core-seconds: OpenFOAM's `ExecutionTime` times the number of ranks.
  J4 times whole-node runs, with `--exclusive`.
- **Disk.** Case directories are deleted after each task. Only the force time series (J2)
  and the observables are kept, so the total is under 1 GB.
- **Adapting to another cluster.** Everything site-specific is in `config.sh` (paths,
  partitions, the OpenFOAM container or module, the GROMACS prefix) and in the `sbatch`
  lines of `common/submit_array.sh` and the `submit.sh` of each job (partition names,
  time limits, memory).

## Revision 4 changes (patch v6)

- **J1.** The function object of the `fo_na*` tasks never wrote: it had `writeControl writeTime` with `writeInterval 20000`, i.e. every 20,000th field write. It now writes every iteration. Every grid run also carries the function object (`cd_fo`, `cl_fo`) and records why each linear solve stopped (`info.stop`); 60 new `rt_*` tasks run relTol 0 and 0.01 at n_arc 32; the `fo_na*` tasks record the wall distance of the centroid and of the vertex mean (`wall`, `wall_geometry.json`).
- **J2.** The spin-up time steps gave maximum Courant numbers of 1.3–1.7 on this mesh family (the article's 2D-2 driver diverged at 1.3); they are now 0.65–0.87, a failed spin-up is retried at half the time step (twice), one failed level no longer stops the others, and the largest Courant number of every spin-up and run is recorded. The analysis saves the case specification that `paper/cases.py` reads as case `2D-2h`, replays every target (marking those not admitted), and prints error against the time step, the relaxation time and the window, the Courant numbers and the accounting of the settings. The cost estimate is 300–600 core-hours: about 1.4 µs per cell and outer iteration with a few outer iterations per step gives roughly 60–100 for the three reference runs and 150–250 for the 180 base runs; runs at large Courant numbers that reach the limit of 50 outer iterations cost more, and the 32-core spin-up job holds its cores for about three hours.
- **Earlier results.** The task farm skips any task with a `result.json` and the J2 spin-up skips any mesh with `SPUN_UP`. The J1 and J2 submit scripts therefore refuse to submit over results written before revision 4, and their smoke tests start from empty smoke folders. The J2 spin-up keeps its attempt number in `SPIN_ATTEMPT`, so a spin job resubmitted after a time-out resumes at the time step it was using, from the checkpoint before the newest one. Smoke analyses write to `results_smoke/`.
- **J3.** The analysis loads the QE case directly (`cases.all_cases()` needs SymPy).
- **J4.** Default configuration v1.2.1.
- **setup_env.sh** installs SymPy and pytest, also into an existing environment.
- **collect_job.sh** runs `collect.sh` as a batch job; `collect.sh` also packs `residuals.npz`, `wall_geometry.json`, `SPUN_UP` and `SPIN_FAILED`.
