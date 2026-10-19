#!/bin/bash
# J4: parallel execution, 1-64 ranks (about 15 node-hours; whole nodes for clean timings).
#   bash submit.sh | bash submit.sh smoke
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
SMOKE=""
if [ "${1:-}" = "smoke" ]; then
    sbatch $SB_EXTRA -p test -N 1 -n 2 -t 00:15:00 --mem=8G -J j4-smoke -o logs/smoke_%j.out <<EOF
#!/bin/bash
source "$KNOBS_ROOT/hamilton/config.sh"; activate_python; cd "$HERE"
export KNOBS_SMOKE=1 WORK="$WORK/smoke"
which \$KNOBS_LMP mpirun; ls \$GMX_DOUBLE_PREFIX/bin/gmx_d
python -m farm parallel_hpc one lmp_R2 && python -m farm parallel_hpc one gmx_R2
python -c "import pickle,glob; [print(f, list(pickle.load(open(f,'rb')).values())[:1]) for f in glob.glob('\$WORK/parallel_hpc/*/tensor.pkl')]"
EOF
    exit 0
fi
# one whole node per task, all tasks in one array
N=$(python -m farm parallel_hpc count)
sbatch $SB_EXTRA -p "$SB_WHOLE" -N 1 --exclusive -t 12:00:00 -J j4 --array=0-$((N-1)) \
    -o logs/j4_%a.out <<EOF
#!/bin/bash
source "$KNOBS_ROOT/hamilton/config.sh"; activate_python; cd "$HERE"
python -m farm parallel_hpc array \$SLURM_ARRAY_TASK_ID $N
EOF
echo "when done: python analyse_j4.py | tee ../results/j4.txt"
