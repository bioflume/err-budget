# Heavy-tier studies on Hamilton

These five jobs produce the results of the revision that need more than the
2-core machine used for the rest of the article. Each job:

- sits in its own directory;
- runs as SLURM arrays of independent, restartable tasks;
- writes one small result per task under `$WORK`;
- has an analysis script that turns those results into the numbers and tables of the article.

| Job | What it answers | Referee point | Core-hours | Wall time |
|---|---|---|---|---|
| `J1_dfg2d1` | 2D-1: forces from OpenFOAM's own function object vs. the article's post-processing; drag convergence to n_arc 192; the wall distance of the two cell centres; the whole 840-setting grid with a current OpenFOAM release, with the function object and the stopping reason of every linear solve; relTol 0 and 0.01 at n_arc 32 | 7, 26 (rev. 4: relTol, centroid offset) | ~20 | < 1 day |
| `J2_dfg2d2` | 2D-2 refined to n_arc 96 (reference at 128 and 192), with the averaging window and the relaxation time as parameters, forces from the function object | 8 (rev. 4: vertex-mean bias) | 300–600 | 1–3 days |
| `J3_qe` | the Quantum ESPRESSO grid with a current release (no `snprintf` shim) | 26 | ~3 | hours |
| `J4_parallel` | LAMMPS at 262,144 ions and GROMACS at 13,824 waters on 1–64 ranks: does communication change the cost ordering or the recommendation? | 11 | ~15 node-hours | < 1 day |
| `J5_abc3d` | **held-out 3D test**: decaying ABC flow with an exact solution, run with the frozen program as pre-registered in `../PREREGISTRATION.md` | 3 | ~60 | 1 day |

These numbers are estimates. They come from per-cell timings on the article's
machine, about 1.4 µs per cell per outer iteration for pimpleFoam. Hamilton's
cores are faster.

## 1. One-off setup

On a login node, about 30–60 minutes:

```bash
tar xzf knobs-1.2.1.tar.gz && cd knobs-1.2.1/hamilton
nano config.sh          # check WORK (scratch), SB_EXTRA (account, if needed), OF_MODE
bash setup_env.sh
```

`setup_env.sh` does four things:

1. Installs micromamba in `~/bin`.
2. Creates the conda environment `knobs`: Python, NumPy, SciPy, a current Quantum ESPRESSO, LAMMPS with OpenMPI, scikit-fem, pyamg and the knobs package.
3. Pulls the OpenCFD container `opencfd/openfoam-default:2406` into `$WORK`.
4. Builds GROMACS 2024 in double precision.

The conda-forge GROMACS builds are single precision, and their round-off floor
(Section 8.3) is too high for the targets of J4.

**OpenFOAM: container or module.** To use Hamilton's own module instead of the
container, set `OF_MODE=module` and `OF_MODULE=` to the name that
`module avail openfoam` shows. The drivers then run in the module's environment.

## 2. Smoke tests (15 minutes each, test partition)

```bash
for J in J1_dfg2d1 J2_dfg2d2 J3_qe J4_parallel J5_abc3d; do (cd $J && bash submit.sh smoke); done
# when they have run:
grep -h "SMOKE TEST" */logs/smoke_*.out
```

Each smoke test runs one or two tiny tasks through the whole pipeline. The
important checks are:

- **J1 smoke test:** that OpenFOAM's `forceCoeffs` output is read. If it fails with "no function-object output", send me `J1_dfg2d1/logs/smoke_*.out`. The parser lives in `common/ofmodern.py` (`read_coefficients`).
- **J2 smoke test:** that `mapFields`, `decomposePar` and a parallel run work in the container.

## 3. Submit

```bash
bash run_all.sh            # submits J1-J5; prints the job ids
squeue -u $USER            # progress
bash status.sh             # tasks done / with errors, per job
```

Every task writes `result.json` when it finishes. Any array can simply be
resubmitted (`bash submit.sh` in its directory): finished tasks are skipped,
and tasks that failed (`error.txt` beside them) are retried.

J2 submits its spin-ups first (one 32-core job, coarse to fine), then the runs
in groups by core count, each group depending on the spin-ups.

## 4. Analyse and send back

```bash
bash collect_job.sh        # runs collect.sh as a batch job (the J1 and J2 analyses take 1-2 hours)
```

This runs the five analysis scripts and writes `results/j1.txt` … `results/j5.txt`
and the tensors (`results/*.pkl`, `results/*.json`). It then packs them into
`hamilton_results.tar.gz`, a few MB. Copy that file to the machine with the article's
repository and run

```bash
python paper/exp60_hamilton.py hamilton_results.tar.gz
```

which unpacks it, replays the held-out problem H3 with the frozen program (from the git
tag, or from `paper/frozen_src` in a release) and with version 1.2.1, regenerates the
tables, and prints the numbers that replace the `[J1]`…`[J5]` placeholders of the article
and the response letter.

## Notes

- **Held-out test (J5).** J5 is the held-out test of the frozen program (version 1.2) and, by the amendment of 11 October 2026, of version 1.2.1. Its problem, parameters, targets and grading rule are fixed in `../PREREGISTRATION.md`, together with the SHA-256 of the frozen source. Do not change `abc_hpc.py` or the targets before running it.
- **Costs.** Costs are core-seconds: OpenFOAM's `ExecutionTime` times the number of ranks. J4 times whole-node runs, with `--exclusive`.
- **Disk.** Case directories are deleted after each task. Only the force time series (J2) and the observables are kept, so the total is under 1 GB.

## Revision 4 changes (patch v6)

- **J1.** The function object of the `fo_na*` tasks never wrote: it had `writeControl writeTime` with `writeInterval 20000`, i.e. every 20,000th field write. It now writes every iteration. Every grid run also carries the function object (`cd_fo`, `cl_fo`) and records why each linear solve stopped (`info.stop`); 60 new `rt_*` tasks run relTol 0 and 0.01 at n_arc 32; the `fo_na*` tasks record the wall distance of the centroid and of the vertex mean (`wall`, `wall_geometry.json`).
- **J2.** The spin-up time steps gave maximum Courant numbers of 1.3–1.7 on this mesh family (the article's 2D-2 driver diverged at 1.3); they are now 0.65–0.87, a failed spin-up is retried at half the time step (twice), one failed level no longer stops the others, and the largest Courant number of every spin-up and run is recorded. The analysis saves the case specification that `paper/cases.py` reads as case `2D-2h`, replays every target (marking those not admitted), and prints error against the time step, the relaxation time and the window, the Courant numbers and the accounting of the settings. The cost estimate is 300–600 core-hours: about 1.4 µs per cell and outer iteration with a few outer iterations per step gives roughly 60–100 for the three reference runs and 150–250 for the 180 base runs; runs at large Courant numbers that reach the limit of 50 outer iterations cost more, and the 32-core spin-up job holds its cores for about three hours.
- **Earlier results.** The task farm skips any task with a `result.json` and the J2 spin-up skips any mesh with `SPUN_UP`. The J1 and J2 submit scripts therefore refuse to submit over results written before revision 4 (move `$WORK/dfg2d1_hpc`, `$WORK/dfg2d2_hpc` and `$WORK/dfg2d2_cases` aside), and their smoke tests start from empty smoke folders. The J2 spin-up keeps its attempt number in `SPIN_ATTEMPT`, so a spin job resubmitted after a time-out resumes at the time step it was using, from the checkpoint before the newest one. Smoke analyses write to `results_smoke/`.
- **J3.** The analysis loads the QE case directly (`cases.all_cases()` needs SymPy).
- **J4.** Default configuration v1.2.1.
- **setup_env.sh** installs SymPy and pytest, also into an existing environment.
- **collect_job.sh** runs `collect.sh` as a batch job; `collect.sh` also packs `residuals.npz`, `wall_geometry.json`, `SPUN_UP` and `SPIN_FAILED`.
