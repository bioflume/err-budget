#!/bin/bash
# J3: Quantum ESPRESSO with a current release (about 3 core-hours).
#   bash submit.sh | bash submit.sh smoke
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
if [ "${1:-}" = "smoke" ]; then
    sbatch $SB_EXTRA -p test -n 1 -t 00:15:00 --mem=4G -J j3-smoke -o logs/smoke_%j.out <<EOF
#!/bin/bash
source "$KNOBS_ROOT/hamilton/config.sh"; activate_python; cd "$HERE"
export KNOBS_SMOKE=1 WORK="$WORK/smoke"
pw.x -h 2>&1 | head -2
python -m farm qe_hpc array 0 1 && python -m farm qe_hpc status
grep -h '"energy"' \$WORK/qe_hpc/g_*/result.json | head -2 && echo "SMOKE TEST PASSED" || echo "SMOKE TEST FAILED"
EOF
    exit 0
fi
bash ../common/submit_array.sh qe_hpc "$HERE" j3 1 12:00:00 4G 40 ""
echo "when done: python analyse_j3.py | tee ../results/j3.txt"
