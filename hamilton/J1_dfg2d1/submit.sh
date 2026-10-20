#!/bin/bash
# J1: DFG 2D-1 with a current OpenFOAM release (about 20 core-hours).
#   bash submit.sh          # all tasks (840 grid runs, 10 function-object levels, 60 relTol runs)
#   bash submit.sh smoke    # 15-minute check on the test partition
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
if [ "${1:-}" = "smoke" ]; then
    sbatch $SB_EXTRA -p test -n 1 -t 00:15:00 --mem=4G -J j1-smoke -o logs/smoke_%j.out <<EOF
#!/bin/bash
source "$ERRBUDGET_ROOT/hamilton/config.sh"; activate_python; load_openfoam; cd "$HERE"
export ERRBUDGET_SMOKE=1 WORK="$WORK/smoke"
rm -rf "\$WORK/dfg2d1_hpc"          # results of an earlier smoke test would be skipped, not re-run
python -m farm dfg2d1_hpc array 0 1
python -m farm dfg2d1_hpc status
python analyse_j1.py
R=\$WORK/dfg2d1_hpc
ok=1
grep -q '"cd_fo": [-0-9]' \$R/fo_na8/result.json || { echo "no function-object output in fo_na8"; ok=0; }
grep -q '"cd_fo": [-0-9]' \$R/grid_0_1_2_0/result.json || { echo "no function-object output in a grid run"; ok=0; }
grep -q '"rel": [0-9]' \$R/grid_0_1_2_0/result.json || { echo "no stopping-criterion counts in a grid run"; ok=0; }
ls \$R/rt_na8_rc0.001_t1e-06_rel0/residuals.npz >/dev/null || { echo "no residual history in the relTol run"; ok=0; }
cat \$R/*/error.txt 2>/dev/null | tail -20
[ \$ok = 1 ] && echo "SMOKE TEST PASSED" || echo "SMOKE TEST FAILED"
EOF
    exit 0
fi
# results written by the driver before revision 4 (no function object, no stopping counts)
# would be kept by the task farm: refuse to submit over them
python - <<'PYEOF' || exit 1
import glob, json, os, sys
R = os.path.join(os.environ["WORK"], "dfg2d1_hpc")
stale = []
for f in glob.glob(os.path.join(R, "*", "result.json")):
    tid, r = os.path.basename(os.path.dirname(f)), json.load(open(f))
    if (tid.startswith("grid_") and "stop" not in r.get("info", {})) or (tid.startswith("fo_na") and "wall" not in r):
        stale.append(f)
if stale:
    sys.exit(f"{len(stale)} results in {R} come from the driver before revision 4 (e.g. {stale[0]}).\n"
             f"Move the folder aside (mv {R} {R}.old) and submit again.")
PYEOF
S=../common/submit_array.sh
N1=$(python -m farm dfg2d1_hpc count nprocs=1)
bash $S dfg2d1_hpc "$HERE" j1-np1 1 12:00:00 4G 60 "" nprocs=1
bash $S dfg2d1_hpc "$HERE" j1-np8 8 2-00:00:00 32G 3 "" nprocs=8
echo "submitted $N1 one-core tasks in 60 array elements, and 3 eight-core tasks"
echo "when done: python analyse_j1.py | tee ../results/j1.txt"
