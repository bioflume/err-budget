#!/bin/bash
# J1: DFG 2D-1 with a current OpenFOAM release (about 15 core-hours).
#   bash submit.sh          # all tasks
#   bash submit.sh smoke    # 15-minute check on the test partition
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
if [ "${1:-}" = "smoke" ]; then
    sbatch $SB_EXTRA -p test -n 1 -t 00:15:00 --mem=4G -J j1-smoke -o logs/smoke_%j.out <<EOF
#!/bin/bash
source "$KNOBS_ROOT/hamilton/config.sh"; activate_python; load_openfoam; cd "$HERE"
export KNOBS_SMOKE=1 WORK="$WORK/smoke"
python -m farm dfg2d1_hpc array 0 1
python analyse_j1.py | head -8
grep -q '"cd_fo"' \$WORK/dfg2d1_hpc/fo_na8/result.json && echo "SMOKE TEST PASSED (function object read)" || echo "SMOKE TEST FAILED: no function-object output"
EOF
    exit 0
fi
S=../common/submit_array.sh
N1=$(python -m farm dfg2d1_hpc count nprocs=1)
bash $S dfg2d1_hpc "$HERE" j1-np1 1 12:00:00 4G 60 "" nprocs=1
bash $S dfg2d1_hpc "$HERE" j1-np8 8 2-00:00:00 32G 3 "" nprocs=8
echo "submitted $N1 one-core tasks in 60 array elements, and 3 eight-core tasks"
echo "when done: python analyse_j1.py | tee ../results/j1.txt"
