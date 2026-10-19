# DFG benchmark 2D-2 (Re = 100) with OpenFOAM: running the build

The build is 243 independent OpenFOAM runs:

- the 240-setting grid: 5 meshes × 4 time steps × 4 outer tolerances × 3 pressure-solver tolerances;
- 3 reference runs on a finer mesh and a smaller time step.

Each run first needs its mesh level spun up once (20 s of flow from rest). Every
finished run is written to its own JSON file in `work_dfg2d2/results/`. Anything
already done is skipped, so you can stop and resubmit at any point.

Spin-ups and runs also save checkpoints as they go: a spin-up every second of flow, a
run every 0.5 s while relaxing and at every sample while measuring. A job killed by a
time limit or a node failure resumes from its last checkpoint when it is resubmitted.
The CPU time it reports is read from the solver logs up to that checkpoint, so it is
that of an uninterrupted run.

## What it costs

Measured on one core of a 2-vCPU cloud VM. A cluster core is usually faster.

| part | CPU time |
|---|---|
| spin-ups (6 meshes, 20 s of flow each) | about 10 h in total; n_arc 48 and 64 take about 3–4 h each |
| grid (240 runs, 4 s of flow each) | about 25 h in total; the longest single runs (n_arc 48, dt 5e-4) take 1–2 h |
| references (3 runs) | about 8 h in total; the longest (n_arc 64, dt 2.5e-4) takes about 3 h |

That is about 40 CPU-hours in all, measured on this VM with three processes sharing two
cores. With 40 array tasks at a time on a cluster it finishes in about 6 hours of wall
time, most of it the n_arc 64 spin-up followed by the reference run.
Memory is below 1 GB per task. The disk peak is about 1 GB per running task, because
fields are written at 500 sampling times and deleted after the forces are extracted.

## Option A: SLURM cluster

1. Copy the package to the cluster:
   ```
   tar xzf knobs-1.1.1.tar.gz
   ```
2. Make Python ready: create a Python environment with numpy and scipy, then run
   `pip install -e knobs-1.1.1`.
3. Make OpenFOAM available. Any ESI-OpenCFD version from v1912 on works: `module load …`, then `source $FOAM_BASH`.
4. Edit `paper/dfg2d2_job/env.sh`. It holds the module and venv lines, plus `KNOBS_DFG2_WORK` if you
   want the work directory on scratch.
5. Add your site's account and partition lines to `slurm_spin.sh` and `slurm_tasks.sh`.
   For example, `#SBATCH --account=...` and `#SBATCH -p ...`.
6. Submit:
   ```
   cd knobs-1.1.1/paper/dfg2d2_job
   ./submit.sh
   ```
   This submits the spin-up array (6 jobs). The task array (243 jobs, at most 40 at once)
   then starts after all spin-ups have succeeded.
7. Monitor progress:
   ```
   ./monitor.sh
   ```
   It shows tasks done, meshes spun up, the estimated worker-hours left, the queue, and logs that contain errors.
8. If tasks time out or fail, run `sbatch slurm_tasks.sh` again from the `dfg2d2_job` directory
   (the scripts find the package through the directory they were submitted from). Finished tasks are skipped.
9. When `monitor.sh` reports 243/243, run:
   ```
   ./finish.sh
   ```
   This writes `paper/data/dfg2d2_tensor.pkl` and `paper/figures/sec8x_dfg2d2.txt`.
10. Send me those two files, and I will add the results to the paper.

## Option B: one machine

```
cd knobs-1.1.1/paper/dfg2d2_job
./run_local.sh 8      # 8 parallel workers (default: all cores)
./monitor.sh
./finish.sh           # when done
```

Workers coordinate through file locks, so starting more workers later is safe.
If the machine reboots, run `./run_local.sh` again.

## Notes

- **Spin-up settings.** Spin-ups use dt = 1e-3 up to n_arc 32 and 5e-4 on the finer
  meshes, which keeps the maximum Courant number below one: at 1.3 the n_arc 48 spin-up
  diverged. The n_arc 64 spin-up runs in PISO mode (one outer iteration), because the
  impulsive start from rest makes the unrelaxed outer iteration diverge in its first
  time step there.
- **Failed runs are expected.** About 30 grid settings with large time steps on the
  finer meshes diverge (maximum Courant number above about 1.5). The analysis treats them
  as failed runs, as the program does.

- **Why OpenFOAM v1912 has no function objects here.** The Ubuntu package of OpenFOAM v1912 cannot construct function objects. So the
  forces are computed from written fields (`openfoam_dfg2d2.py`), the same way on
  every installation. On a cluster OpenFOAM, the numbers therefore come from the same
  code path. Only the CPU times will differ from those of the 2D-1 study.
- **Why the reference comes from OpenFOAM itself.** The published reference values for 2D-2 are intervals: c_D,max in
  [3.22, 3.24], c_L,max in [0.99, 1.01], St in [0.295, 0.305]. That is too coarse to
  grade targets of 1e-3. The analysis therefore uses the finest OpenFOAM run as the
  reference, measures its uncertainty from the two other reference runs, and admits
  only targets at least ten times that uncertainty.
