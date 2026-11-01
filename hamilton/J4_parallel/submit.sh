#!/bin/bash
# J4: parallel execution, 1-64 ranks (about 15 node-hours; whole nodes for clean timings).
#   bash submit.sh | bash submit.sh smoke
#
# November 2026: every run of October failed within a second and was stored as a finished
# task with no value; such results are refused (move them aside first).  Each task asks for
# 64 tasks on its node, so that mpirun has a slot for every rank, and LAMMPS's own
# messages are kept when a run fails (paper/exp55_parallel.py, LAST_ERROR).
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
SMOKE=""
if [ "${1:-}" = "smoke" ]; then
    sbatch $SB_EXTRA -p test -N 1 -n 2 -t 00:15:00 --mem=8G -J j4-smoke -o logs/smoke_%j.out <<EOF
#!/bin/bash
source "$ERRBUDGET_ROOT/hamilton/config.sh"; activate_python; cd "$HERE"
export ERRBUDGET_SMOKE=1 WORK="$WORK/smoke"
rm -rf "\$WORK/parallel_hpc"     # an earlier smoke test would be skipped, not re-run
which \$ERRBUDGET_LMP mpirun; ls \$GMX_DOUBLE_PREFIX/bin/gmx_d
python -m farm parallel_hpc one lmp_R2; python -m farm parallel_hpc one gmx_R2
cat \$WORK/parallel_hpc/*_R2/error.txt 2>/dev/null | tail -n 20
python -c "import pickle,glob; [print(f, list(pickle.load(open(f,'rb')).values())[:1]) for f in glob.glob('\$WORK/parallel_hpc/*/tensor.pkl')]"
python -c "import glob,math,pickle,sys; f=sorted(glob.glob('\$WORK/parallel_hpc/*_R2/tensor.pkl')); sys.exit(0 if len(f) == 2 and all(any(math.isfinite(v[1]) for v in pickle.load(open(x,'rb')).values()) for x in f) else 1)" && echo "SMOKE TEST PASSED" || echo "SMOKE TEST FAILED"
EOF
    exit 0
fi

# refuse to submit twice, or over the failed results of October
if squeue -h -u "$USER" -o "%j" 2>/dev/null | grep -qx 'j4-grid'; then
    echo "J4 jobs are already queued or running (squeue -u $USER); submitting again would run tasks twice."
    exit 1
fi
python - <<'PYEOF' || exit 1
import glob, json, math, os, pickle, sys
W = os.environ["WORK"]
bad = []
for f in glob.glob(os.path.join(W, "parallel_hpc", "*", "tensor.pkl")):
    if not any(math.isfinite(v[1]) for v in pickle.load(open(f, "rb")).values()):
        bad.append(os.path.dirname(f))
for f in glob.glob(os.path.join(W, "parallel_hpc", "*", "ref.json")):
    if not math.isfinite(json.load(open(f)).get("QS", float("nan"))):
        bad.append(os.path.dirname(f))
if bad:
    sys.exit(f"{len(bad)} J4 results in {W}/parallel_hpc have no value (failed runs):\n  "
             + "\n  ".join(sorted(bad)) + "\nMove the folder aside (do not delete it), e.g.\n"
             f"  mkdir -p {W}/before_rerun && mv {W}/parallel_hpc {W}/before_rerun/\n"
             "and submit again.")
PYEOF
if [ "$(python -m farm parallel_hpc left)" = "0" ]; then
    echo "J4: all tasks done"
    exit 0
fi
[ -f "$WORK/J4_SUBMITTED" ] || date +%F > "$WORK/J4_SUBMITTED"

# one whole node per task, all tasks in one array (finished tasks exit at once)
N=$(python -m farm parallel_hpc count)
sbatch $SB_EXTRA -p "$SB_WHOLE" -N 1 --exclusive -n 64 -t 12:00:00 -J j4-grid --array=0-$((N-1)) \
    -o logs/j4_%a.out <<EOF
#!/bin/bash
source "$ERRBUDGET_ROOT/hamilton/config.sh"; activate_python; cd "$HERE"
python -m farm parallel_hpc array \$SLURM_ARRAY_TASK_ID $N
EOF
echo "monitor: squeue -u \$USER;  bash ../status.sh  (node-hours so far: bash ../coreh.sh j4)"
