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

## November 2026: finishing J2 and J4 (read this first)

The first campaign (20–31 October 2026) finished J1 and J5 (H3). J2, J3 and J4 have to be
redone; J3 now runs on the article's machine, so only J2 and J4 come back to Hamilton.
**J1 and J5 are final: do not submit them again, and leave `$WORK/dfg2d2_cases/mesh_24`,
`$WORK/abc_hpc` and `$WORK/dfg2d1_hpc` as they are.**

What went wrong, and what this version changes:

- **Container launches.** Open MPI writes its session directory under `$TMPDIR`, which
  SLURM on Hamilton sets to node-local `/local/...`; `/local` does not exist in the
  OpenFOAM container, so every parallel run failed with "mkdir /local: Read-only file
  system" until it was bound in on 22 October. `config.sh` now binds `/local` (and uses
  `/tmp` in the container where `/local` is absent).
- **Paths of codes in their own environments** (`config.local.sh`, written in October for
  LAMMPS and Quantum ESPRESSO): `config.sh` reads it before its defaults, and accepts the
  old names (`KNOBS_LMP`, `KNOBS_PW`, `KNOBS_MPIRUN`).
- **J2, the inlet.** `mapFields -consistent` maps the values of every patch, the
  fixedValue inlet included. Every mesh that was mapped from a coarser one (n_arc 32 and
  finer, the references included) therefore ran with the coarser mesh's inlet profile,
  interpolated, instead of the analytic one: the flow rate differs by up to 0.4 % (n_arc 48
  mapped from 32), and the n_arc 48 runs differ from those of revision 3 by 4.7e-3 in the
  Strouhal number. `spin()` now writes the analytic inlet back after mapping, checks it and
  records it in `mesh_N/INLET_ANALYTIC`; the driver, `submit.sh` and `analyse_j2.py` refuse
  meshes and results from before the fix. Only n_arc 24 (spun up from rest) is unaffected.
- **J2, the spin-ups.** One job per mesh, on the cores its runs use, with a one-day limit,
  each after the coarser one; the runs of a mesh start when its spin-up has succeeded, and
  SLURM cancels whatever depends on a failed job. A spin-up that fails before its first time
  step (decomposePar, mpirun) no longer halves the time step: in October such failures used up
  both halvings, and the meshes were spun up at a quarter of the registered step. The 32-core
  spin-up job that held its cores for seven days (5,376 of J2's 6,891 core-hours) cannot recur.
- **J2, slowed-down runs.** The spin-up of n_arc 192 on 32 ranks ran at about 4 s per
  linear-solver iteration, against 2–5 ms for the coarser meshes in the same job, and the
  reference at n_arc 128 and dt 2.5e-4 on 16 ranks took about 1,250 s per time step, against
  0.9 s for the one at dt 1.25e-4 on another node; both used up their time limits. mpirun now
  runs with `--bind-to none`, and J2 stops a run whose linear solves and their iterations
  take more than 0.1 s each over 15 minutes (normal: 2–6 ms; `ERRBUDGET_OF_MAX_S_PER_ITER`,
  `ERRBUDGET_OF_WATCH_WINDOW`). Such a stop is an error, not a failed setting: the task is
  retried at the next submission and resumes from its last checkpoint, and the work it
  redoes is not charged to its cost. A resumed run that fails before its first time step
  (a launch failure, or a checkpoint cut off mid-write) is an error too, retried from the
  checkpoint before.
- **J4.** LAMMPS works with the environment of `config.local.sh` (checked by the probe of
  1 November). GROMACS needs the C++ runtime of the compiler that built it: every job now
  loads `gcc/11.2` (`ERRBUDGET_MODULES`). Each task asks for 64 tasks on its whole node, so
  that mpirun has a slot for every rank; every LAMMPS task first runs one cheap setting and
  stops with LAMMPS's own messages if it fails; a LAMMPS run on fewer processes than ranks (an
  mpirun that does not match the MPI of `lmp`) is a failed run; `submit.sh` refuses the
  failed results of October.
- **Task farm.** An `error.txt` from an earlier attempt is kept as `error_prev.txt`, so a
  time-out is no longer reported with an old error; every attempt is logged in
  `attempts.txt`; `status.sh` counts tasks started and unfinished and prints the core-hours
  of the reruns (`coreh.sh`).
- **J3 and J4.** A run that produced no output at all (pw.x, LAMMPS or GROMACS did not
  start) is now an error, retried, with the program's own error message, instead of a
  stored NaN; the smoke tests can fail (the J3 one could not).

### Steps

1. **Code.** Use this repository at the commit of these changes (`git log -1`), not a copy
   of an older version. In `hamilton/config.sh` set `WORK=/nobackup/$USER/knobs` and
   `ERRBUDGET_ENV=knobs` (the environment and scratch folder of October), copy
   `config.local.sh` from the October copy into `hamilton/`, and carry over any other
   setting you changed in October (OpenFOAM mode or module, decomposition method, account).
2. **Move the old J2 and J4 state aside** (never delete it):
   ```bash
   W=/nobackup/$USER/knobs; mkdir -p $W/before_inlet_fix $W/before_rerun
   mv $W/dfg2d2_hpc $W/before_inlet_fix/
   for m in 32 48 64 96 128 192; do mv $W/dfg2d2_cases/mesh_$m $W/before_inlet_fix/; done
   mv $W/dfg2d2_cases/SPIN_FAILED $W/before_inlet_fix/ 2>/dev/null
   mv $W/parallel_hpc $W/before_rerun/
   ```
3. **Smoke tests** (15 minutes each; the J2 one maps its finer mesh and runs it on two
   ranks, as in production): `bash submit.sh smoke` in `J2_dfg2d2` and in `J4_parallel`;
   then `grep "SMOKE TEST" */logs/smoke_*.out` must say PASSED for both.
4. **Submit**: `bash submit.sh` in each (they refuse if their jobs are already queued). J2
   prints one spin-up job per mesh (32 to 192) and one array per mesh; J4 one array of whole
   nodes.
5. **Monitor** daily with `bash status.sh`. Each spin-up log (`J2_dfg2d2/logs/spin<N>_*.out`)
   should say "analytic inlet restored" and end with "spun up". The guard is 2,200 core-hours
   for the J2 rerun: if `status.sh` prints GUARD PASSED, stop and report. A spin-up or a run
   stopped as too slow ("was stopped" in its spin-up log or in an `error.txt`) and time-outs:
   run `bash submit.sh` again once the queue holds no job of that study; finished tasks and
   meshes are skipped, the others resume. If a spin-up fails for another reason, its
   dependants are cancelled: report its log. A task or spin-up that fails twice for the same
   reason: report it.
6. **Collect** only the reruns: `JOBS="2 4" bash collect_job.sh` (Section 5 below).

Expected cost: J2 about 1,500–2,000 core-hours (base runs about 450, references about
600–1,000, spin-ups about 150), 1–3 days; J4 15 whole-node tasks of up to 12 hours.

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

1. **The held-out test (J5) stays as registered.** Do not change `src/errbudget/`,
   `../PREREGISTRATION.md`, `J5_abc3d/abc_hpc.py` or the H3 targets. J5 runs once. Its
   results are reported as they come out, good or bad. The program was called knobs when
   J5 was registered; `abc_hpc.py` and the frozen program keep that name, which is why
   the J5 smoke test also sets `KNOBS_SMOKE`.
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
git clone <this repository> errbudget-paper && cd errbudget-paper/hamilton
nano config.sh        # WORK (scratch, default /nobackup/$USER/errbudget); SB_EXTRA (account, if
                      # yours needs one); OF_MODE (container or module, see below)
bash setup_env.sh     # safe to re-run: finished steps are skipped
```

`setup_env.sh` installs micromamba in `~/bin`; creates the conda environment `errbudget`
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
python exp60_hamilton.py /path/to/hamilton_results.tar.gz
bash replay_all.sh
```

`exp60_hamilton.py` unpacks the archive into `hamilton/results`, prints J1–J5, replays
H1–H3 with the frozen program and with the released one (the line `frozen program
verified` must appear) and regenerates the tables. `replay_all.sh` then replays every
result, the refined 2D-2 study (case `2D-2h`) included. If `hamilton/results/j2_tensor.pkl`
lacks the case specification, `pooled.py` skips case 2D-2h and says so: re-run the patched
`J2_dfg2d2/analyse_j2.py` on Hamilton. The small files of `hamilton/results` (everything but
`raw.tar.gz`, 1 MB) are in the repository, so that the replay includes the refined 2D-2
study from a clone; the raw outputs (`raw.tar.gz`, 47 MB) go to the Zenodo archive
(`../DATA.md`).

## 7. What to send back

1. `hamilton_results.tar.gz` (from Section 5) and, after Section 6, `errbudget_v6_results.tar.gz`.
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
