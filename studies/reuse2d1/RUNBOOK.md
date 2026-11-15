# Runbook: the reuse study on Hamilton

Lines marked **(Mac)** run in a Terminal on the Mac; everything else in a Hamilton session
(`ssh grbj73@hamilton8.dur.ac.uk`). It uses the airfoil study's OpenFOAM module
(`openfoam/v2606`) and Python environment (`~/venvs/airfoil`); nothing new is installed.

| Step | What | Partition | Core-hours | Wall time |
|---|---|---|---|---|
| check | 3 grid runs end to end | `test` | < 0.1 | 15 min |
| grid | 840 settings × 13 Reynolds numbers, 91 array tasks | `shared` | about 30–50 | hours |
| ref | 39 reference runs | `shared` | about 5–10 | hours |

Rules: do not change numerical settings (lists, members, meshes, tolerances, schemes);
what may change is what stops a job from running (module names, partition, account,
memory, time limits). A failed or unconverged run is a result; never delete files in
`runs/`. Stop and report if the check fails, if `status.sh` prints GUARD PASSED, or if a
task is given up after three errors.

## 1. Upload (Mac)

```bash
RE="/Users/grbj73@durham.ac.uk/.local/share/openresearch/worktrees/2ed68104-cb7b-456e-999e-56a5d09d1af9/chat_7c2bdc9c-696d-ead8-f4ba-6e51075f2c4f/reuse2d1"
shasum -a 256 "$RE/dist/reuse2d1-batch.tar.gz"; cat "$RE/dist/SHA256SUMS"     # the two must agree
scp "$RE/dist/reuse2d1-batch.tar.gz" grbj73@hamilton8.dur.ac.uk:/nobackup/grbj73/
```

## 2. Unpack (Hamilton)

```bash
W=/nobackup/$USER/reuse2d1
mkdir -p $W/code && tar xzf /nobackup/$USER/reuse2d1-batch.tar.gz -C $W/code && cat $W/code/VERSION
# the airfoil study's local settings (module name, OpenFOAM bashrc), if it has any
A=/nobackup/$USER/airfoil/code/slurm/config.local.sh
[ -f $A ] && grep -E 'AIRFOIL_(OF_MODULE|OF_BASHRC|PRE_MODULES)=' $A | sed 's/AIRFOIL_/REUSE_/g' > $W/code/slurm/config.local.sh
cat $W/code/slurm/config.local.sh 2>/dev/null
```

## 3. The check job (`test` partition, 15 minutes)

```bash
bash $W/code/slurm/submit.sh check
squeue -u $USER                                    # wait until the job has gone
tail -n 12 $(ls -t $W/logs/check_*.out | head -1)  # last line: CHECK PASSED
```

| CHECK FAILED because | Do |
|---|---|
| `simpleFoam not on PATH`, `cannot load`, `no virtual environment` | the module name or the venv: as for the airfoil study (its RUNBOOK.md, steps 2–3) |
| `version ..., not v2606` | choose the v2606 module |
| anything else | send me the check log and `$W/tasks/*.error-*.txt` |

## 4. Submit

```bash
bash $W/code/slurm/submit.sh grid      # 91 tasks in two arrays (4 h / 2 GB and 12 h / 3 GB)
bash $W/code/slurm/submit.sh ref       # 39 tasks in two arrays (6 h / 4 GB and 12 h / 6 GB)
```

## 5. Watch and resubmit

```bash
bash $W/code/slurm/status.sh           # tasks and runs, the queue, core-hours (guard 200), errors
bash $W/code/slurm/submit.sh grid; bash $W/code/slurm/submit.sh ref   # again after failures or time-outs
```

The same `submit.sh` line resubmits: finished tasks and tasks still running are skipped,
an interrupted task goes on from its last stored run, and a task with an error file is
tried again, at most three times. After an OUT_OF_MEMORY or TIMEOUT that comes back:
`REUSE_MEM=8G REUSE_TIME=1-00:00:00 bash $W/code/slurm/submit.sh grid`.

## 6. Pack and copy back

When `status.sh` shows every task done for grid and ref:

```bash
bash $W/code/slurm/pack.sh
scp grbj73@hamilton8.dur.ac.uk:/nobackup/grbj73/reuse2d1/packs/reuse2d1-<date>.tar.gz "$RE/results/"   # (Mac)
```

The analysis runs on the Mac (`analysis/reuse.py`).
