#!/bin/bash
# J5: held-out 3D test, decaying ABC flow (about 60 core-hours).
#   bash submit.sh | bash submit.sh smoke
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
if [ "${1:-}" = "smoke" ]; then
    sbatch $SB_EXTRA -p test -n 1 -t 00:15:00 --mem=4G -J j5-smoke -o logs/smoke_%j.out <<EOF
#!/bin/bash
source "$ERRBUDGET_ROOT/hamilton/config.sh"; activate_python; load_openfoam; cd "$HERE"
# abc_hpc.py is fixed by the pre-registration and reads the old name KNOBS_SMOKE
export ERRBUDGET_SMOKE=1 KNOBS_SMOKE=1 WORK="$WORK/smoke"
python -m farm abc_hpc array 0 1 && grep -h '"ek"' \$WORK/abc_hpc/*/result.json | head -1 && echo "SMOKE TEST PASSED" || echo "SMOKE TEST FAILED"
EOF
    exit 0
fi
S=../common/submit_array.sh
for spec in "1 12:00:00 4G 36" "2 12:00:00 8G 36" "4 24:00:00 16G 36" "8 2-00:00:00 32G 36"; do
    set -- $spec
    N=$(python -m farm abc_hpc count nprocs=$1)
    [ "$N" = "0" ] && continue
    bash $S abc_hpc "$HERE" "j5-np$1" "$1" "$2" "$3" "$(( N < $4 ? N : $4 ))" "" "nprocs=$1"
done
echo "when done: python analyse_j5.py | tee ../results/j5.txt"
