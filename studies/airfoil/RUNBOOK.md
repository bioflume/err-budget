# Runbook: the airfoil study on Hamilton

The commands are meant to be pasted in order. Lines marked **(Mac)** run in a Terminal on
the Mac; everything else runs in a Hamilton session (`ssh grbj73@hamilton8.dur.ac.uk`).
Each stage ends with one tarball copied back to the Mac. Stop at the end of each stage:
the next one starts after the tarball has been checked.

| Stage | What | Partition | Core-hours (predicted) | Wall time |
|---|---|---|---|---|
| check | 2 tiny runs end to end | `test` | < 0.1 | 15 min |
| 0 | 26 runs: cost, convergence, hardest design, far field | `shared` | 5 (guard 8) | a few hours |
| 1 | 540-run grid + 5 reference runs at the hardest design | `shared` | 80 (guard 400) | about a day |
| 2 | 25 designs × 4 settings | `shared` | 130 (guard 500) | about a day |

Background: PROTOCOL.md (what and why), COST.md (the estimate).

## Rules

1. **Do not change numerical settings** (lists of values, targets, meshes, tolerances,
   schemes, the reference). What may change is what stops a job from running: module names,
   partitions, account, memory and time limits. Record every change and report it.
2. **A failed run is a result.** A run that diverges or does not converge is stored as such;
   never delete or edit a file in `runs/` to make a failure go away. Files named
   `<tag>.error-<k>.json` are infrastructure errors: fix the cause and resubmit.
3. **Stop and report** when the check job fails for a reason other than the module name;
   when `status.sh` prints GUARD PASSED; when a run is given up after three errors; when
   the Stage-0 analysis prints STOP.

## 1. Upload (Mac)

`AF` is the `airfoil/` directory of a checkout of branch `study/airfoil`; the line below
names the worktree where it was prepared.

```bash
AF="/Users/grbj73@durham.ac.uk/.local/share/openresearch/worktrees/2ed68104-cb7b-456e-999e-56a5d09d1af9/chat_7c2bdc9c-696d-ead8-f4ba-6e51075f2c4f/airfoil"
shasum -a 256 "$AF/dist/airfoil-batch.tar.gz"; cat "$AF/dist/SHA256SUMS"     # the two must agree
scp "$AF/dist/airfoil-batch.tar.gz" grbj73@hamilton8.dur.ac.uk:/nobackup/grbj73/
```

## 2. Unpack and choose the OpenFOAM module (Hamilton)

```bash
W=/nobackup/$USER/airfoil
mkdir -p $W/code
[ -d $W/code/config ] && cp -r $W/code/config $W/config.bak-$(date +%Y%m%d-%H%M%S)   # later uploads
tar xzf /nobackup/$USER/airfoil-batch.tar.gz -C $W/code
cat $W/code/VERSION
module avail openfoam 2>&1 | grep -i -E "2606|openfoam"
module avail OpenFOAM 2>&1 | grep -i 2606
```

Take the exact name of the ESI **v2606** module from that list (for example `openfoam/v2606`),
and write it into `config.local.sh`, which later uploads do not overwrite:

```bash
echo 'export AIRFOIL_OF_MODULE=openfoam/v2606' > $W/code/slurm/config.local.sh   # <- the name you saw
module show openfoam/v2606 2>&1 | head -30                                       # same name
```

If `module show` says the module needs others loaded first (a compiler, an MPI), add them:
`echo 'export AIRFOIL_PRE_MODULES="gcc/11.2 openmpi/4.1.1"' >> $W/code/slurm/config.local.sh`.

## 3. One-time setup (login node, a few minutes)

```bash
bash $W/code/slurm/setup.sh 2>&1 | tee $W/setup.log
```

It makes the Python environment `~/venvs/airfoil` (the default Python module is 3.9.9;
NumPy, SciPy, PyYAML from pip), loads OpenFOAM, prints `simpleFoam`'s path and version,
copies the airFoil2D tutorial's `fvSolution`, and ends with `check: 2 runs` and
`stage0: 26 runs`.

Hamilton's `openfoam/v2606` module only sets `OPENFOAM_HOME` and prints "Please now execute:
source $OPENFOAM_HOME/OpenFOAM/etc/bashrc". `env.sh` sources that file itself. If it still
says OpenFOAM did not load, give the path explicitly and run setup again:

```bash
module load openfoam/v2606 && ls $OPENFOAM_HOME/OpenFOAM/etc/bashrc
echo "export AIRFOIL_OF_BASHRC=$OPENFOAM_HOME/OpenFOAM/etc/bashrc" >> $W/code/slurm/config.local.sh
bash $W/code/slurm/setup.sh 2>&1 | tee $W/setup.log
```

## 4. The check job (`test` partition, 15 minutes)

```bash
bash $W/code/slurm/submit.sh check
squeue -u $USER                                   # wait until the job has gone
tail -n 40 $(ls -t $W/logs/check_*.out | head -1) # six lines "MESH PASSED", then the last line: CHECK PASSED
```

The job first writes six meshes (the coarsest and finest of the Stage-1 lists for two
corners, the reference and the n_surf 512 mesh) and runs checkMesh on each; a mesh that
fails a check other than the aspect ratio fails the job (PROTOCOL.md §14, A2). Then it runs
two tiny cases end to end. Results made with an older mesh or case set-up (A3) are moved
to `$W/runs/stale/`.

| CHECK FAILED because | Do |
|---|---|
| `simpleFoam not on PATH`, `cannot load` | the module name or AIRFOIL_PRE_MODULES (step 2) |
| `MESH FAILED` | send me the check log and `$W/meshcheck/*.checkMesh.txt`; do not edit the code |
| an error file in `$W/runs` mentions checkMesh or a keyword | send me the check log and `$W/runs/*.error-1.json`; do not edit the code |
| `version ..., not v2606` | the module is another release: choose the v2606 one |
| `did not echo decayControl` | send me the check log (the free-stream turbulence would decay) |
| `tutorial ...` values differ | send me `$W/tutorial_airFoil2D_fvSolution` |

To run the check again, remove only its two runs and resubmit:

```bash
rm -f $W/runs/m2.0_t10_n64_y2e-05_R50_rc1e-03_gt1e-07.* $W/runs/m4.0_t12_n64_y2e-05_R50_rc1e-03_gt1e-07.*
bash $W/code/slurm/submit.sh check
```

## 4b. The confirmation map (before Stage 0; PROTOCOL.md §14, A4)

Amendment A4 narrowed the design space to camber 0–4 %, thickness 8–12 %, and starts every
run finer than the Stage-1 lists (the reference, the n_surf 512 check) from the same
design's n_surf 91 solution, mapped onto its mesh (`mapFields`). After the check job has
passed with this tarball, confirm convergence where it is hardest: the 45 Stage-1 meshes at
residualControl 10⁻⁶ for m 4 %, t 12 % and m 4 %, t 8 %, and the reference of m 4 %, t 12 %
(91 runs; the convergence map of A3 is superseded and is not run again).

```bash
bash $W/code/slurm/submit.sh confirm             # once; about 40 core-hours, at most about 80
bash $W/code/slurm/status.sh confirm             # repeat until no run is running or pending
source $W/code/slurm/env.sh python && cd $W/code && python3 analysis/convmap.py --sweep confirm
bash $W/code/slurm/pack.sh confirm               # then copy the archive to the Mac
```

The 90 mesh runs finish within about two hours; the reference (846,720 cells, about 10,000
iterations after its first run) takes about 11 hours. Pack when the 90 are done, and again
when the reference has finished. Its log starts with "first run m4.0_t12_n91_...: ...,
mapped onto ...".

Submit once. `submit.sh` leaves out runs that are running: while a run executes it holds
`$W/runs/<tag>.running`, and a second copy of it exits at once ("RUNNING ELSEWHERE", driver
exit 5) instead of sharing its case directory.

## 4c. The settling measurement (12 November; no registered rule changes)

The reference of the confirmation map stalled with L/D in a slow, damped oscillation of about
±2 %. Two runs of the reference mesh of m 4 %, t 12 %, at R = 50,000 and at R = 5,000, run to
40,000 iterations without a residual stop or a stall watchdog, each started as in A4 from its
n_surf 91 solution; their whole force history goes to `$W/runs/history/<tag>.csv`.

```bash
bash $W/code/slurm/submit.sh settle              # once; 2 runs, limit 2.5 days each
bash $W/code/slurm/status.sh settle              # done 2 when both have finished
source $W/code/slurm/env.sh python && cd $W/code && python3 analysis/settle.py
bash $W/code/slurm/pack.sh settle                # then copy the archive to the Mac
```

Expected: about 20 core-hours for R = 50,000 and 16 for R = 5,000 at the reference's measured
speed (2.1 µs per cell and iteration), up to about 43 and 35 on slower nodes; the two run side
by side, so about one day of wall time. If a run reaches its time limit, Slurm warns the job
ten minutes before, the run writes its fields and stops ("INTERRUPTED ... resumes from"), and
`status.sh settle` shows it as started: submit the same line again and it continues.

## 5. Stage 0

```bash
bash $W/code/slurm/submit.sh stage0
bash $W/code/slurm/status.sh check stage0         # repeat until "stage0: 26 runs: done 26"
```

`submit.sh` submits one array per resource class (time limit, memory, cores), only for runs
that have no result yet, and records the job ids in `$W/jobs.txt`.

## 6. Check Stage 0, pack, copy back

```bash
source $W/code/slurm/env.sh python && cd $W/code
mkdir -p $W/results
python3 analysis/stage0.py --write-config 2>&1 | tee $W/results/stage0.txt
bash slurm/pack.sh stage0
```

`stage0.py` prints the accounting, the measured cost model, y⁺, the predicted core-hours of
Stages 1–2, the hardest design, and the far-field and numerics checks, and writes
`config/stage1.json`. If it says Stage 0 is not complete, resubmit (step 8) first. The last
lines of `pack.sh` print the copy line; on the Mac:

```bash
mkdir -p "$AF/results"                                                    # (Mac)
scp grbj73@hamilton8.dur.ac.uk:/nobackup/grbj73/airfoil/packs/airfoil-stage0-<date>.tar.gz "$AF/results/"   # (Mac)
```

**Stop here.** Stage 1 starts after the Stage-0 tarball has been checked.

## 7. Monitoring (every stage)

```bash
squeue -u $USER
bash $W/code/slurm/status.sh                      # runs per sweep, queue, core-hours against the guards
sacct -X -u $USER -S $(date -d '3 days ago' +%F) --format=JobID%22,JobName%14,State,Elapsed,AllocCPUS,ReqMem,MaxRSS
tail -n 20 $W/logs/af-stage1_<job>_<index>.out    # one task's log; the driver's last line says what happened
```

"PD (Priority)" in `squeue` means waiting behind other jobs, not failed; `squeue --start -u
$USER` gives an estimate.

## 8. Resubmitting failed, interrupted or timed-out runs

```bash
source $W/code/slurm/env.sh python && cd $W/code
python3 analysis/accounting.py --runs $W/runs stage0          # or stage1 stage1ref, or stage2
bash slurm/submit.sh stage0                                   # submits only the runs without a result
```

The same `submit.sh` line is the resubmission: finished runs and runs still running are
skipped (a running run holds `$W/runs/<tag>.running`), interrupted runs
(time limit; driver exit 3) resume from the fields they wrote, and runs with an error file
are tried again, at most three times in all. After an OUT_OF_MEMORY or a TIMEOUT that keeps
coming back, raise the limit for that submission only:
`AIRFOIL_MEM=8G AIRFOIL_TIME=1-00:00:00 bash slurm/submit.sh stage1`. A run given up after
three errors is listed by `accounting.py`: report it with its error file.

## 9. Stage 1 (after the Stage-0 check)

```bash
cd $W/code && head -8 config/stage1.json            # the hardest design (the radius factor stays 1, A5)
bash slurm/submit.sh stage1                         # 540 runs, three arrays by mesh size
bash slurm/submit.sh stage1ref                      # the 6 reference runs (A5; time limits up to 3 days)
bash slurm/status.sh stage1 stage1ref               # daily; guard 400 core-hours
bash slurm/submit.sh stage1; bash slurm/submit.sh stage1ref   # again after failures or time-outs
```

When both are complete (`accounting.py ... stage1 stage1ref` says COMPLETE twice):

```bash
bash slurm/pack.sh stage1
scp grbj73@hamilton8.dur.ac.uk:/nobackup/grbj73/airfoil/packs/airfoil-stage1-<date>.tar.gz "$AF/results/"   # (Mac)
```

**Stop here.** The Stage-1 analysis writes `config/stage2_settings.json`, and you receive a
new tarball for Stage 2.

## 10. Optional: errbudget live, through its job file (after Stage 1)

```bash
cd $W/code && bash slurm/submit.sh live            # minutes: every run it asks for is stored
tail -n 30 $W/logs/live_*.out                      # errbudget's report
bash slurm/pack.sh live
scp grbj73@hamilton8.dur.ac.uk:/nobackup/grbj73/airfoil/packs/airfoil-live-<date>.tar.gz "$AF/results/"   # (Mac)
```

It must return the setting that the replay of Stage 1 returned for (a). If the far-field
amendment was triggered, `live_tune.py` refuses to run until the radius values of
`errbudget_job.yaml` are multiplied by ten (the new tarball has them).

## 11. Stage 2 (with the new tarball)

Upload and unpack the new tarball as in steps 1 and 2 (the first `cp` line keeps a copy of
the old configuration; `config.local.sh` is not touched), then:

```bash
cd $W/code && cat VERSION && cat config/stage2_settings.json
bash slurm/submit.sh stage2                         # 125 runs (A5: 25 far-field runs), fewer if some exist
bash slurm/status.sh stage2                         # guard 500 core-hours
bash slurm/submit.sh stage2                         # again after failures or time-outs
bash slurm/pack.sh stage2
scp grbj73@hamilton8.dur.ac.uk:/nobackup/grbj73/airfoil/packs/airfoil-stage2-<date>.tar.gz "$AF/results/"   # (Mac)
```

## 12. Analysis on the Mac (for reference; Claude runs it)

```bash
cd "$AF"                                                                            # (Mac)
mkdir -p results/stage1 && tar xzf results/airfoil-stage1-<date>.tar.gz -C results/stage1
PYTHONDONTWRITEBYTECODE=1 /opt/anaconda3/bin/python3 analysis/stage1.py \
    --runs results/stage1/runs --config results/stage1/packs/meta-stage1-<date>/config \
    --out results/stage1/analysis --write-settings
```

`analysis/stage0.py` and `analysis/stage2.py` take the same arguments; `stage2.py` also
`--stage1 results/stage1/analysis/stage1_summary.json`. `results/` is not tracked by git.

## What is where on Hamilton

```
/nobackup/grbj73/airfoil/
  code/       this tarball: driver.py, af/ (geometry, mesh, case), slurm/, analysis/, config/, engine/
  runs/       the result store: <tag>.json and <tag>.log.txt per run, <tag>.error-<k>.json
  cases/      OpenFOAM case directories while they run (deleted after; kept if interrupted)
  logs/       Slurm logs: check_<job>.out, af-<sweep>_<job>_<index>.out
  results/    analysis output; packs/: the tarballs; jobs.txt: every job id submitted
```

Each run keeps two small files (inode quota); case directories and fields are deleted when a
run finishes. A tarball packed while jobs run is a snapshot, and GNU tar may then exit 1
("file changed as we read it"); that is not an error.
