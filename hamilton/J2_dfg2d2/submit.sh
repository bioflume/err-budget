#!/bin/bash
# J2: DFG 2D-2 refined (about 200 core-hours; 1-2 days of queue and wall time).
#   bash submit.sh          # spin-ups, then the runs in groups by core count
#   bash submit.sh smoke    # a 15-minute check (two tiny meshes)
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
if [ "${1:-}" = "smoke" ]; then
    export KNOBS_SMOKE=1
    sbatch --parsable $SB_EXTRA -p test -n 2 -t 00:15:00 --mem=4G -J j2-smoke \
        -o logs/smoke_%j.out --export=ALL,KNOBS_SMOKE=1 <<EOF
#!/bin/bash
source "$KNOBS_ROOT/hamilton/config.sh"; activate_python; load_openfoam; cd "$HERE"
export KNOBS_SMOKE=1 WORK="$WORK/smoke"
python dfg2d2_hpc.py spin && python -m farm dfg2d2_hpc array 0 1 && python -m farm dfg2d2_hpc status
ls \$WORK/dfg2d2_hpc/*/series.npz && echo "SMOKE TEST PASSED" || echo "SMOKE TEST FAILED"
EOF
    exit 0
fi
# 1. spin-ups, coarse to fine (n_arc 192 uses 32 cores)
SPIN=$(sbatch --parsable $SB_EXTRA -p "$SB_SHARED" -N 1 -n 32 -t 2-00:00:00 --mem=64G -J j2-spin \
    -o logs/spin_%j.out <<EOF
#!/bin/bash
source "$KNOBS_ROOT/hamilton/config.sh"; activate_python; load_openfoam; cd "$HERE"
python dfg2d2_hpc.py spin
EOF
)
echo "spin-ups: job $SPIN"
# 2. runs, grouped by the number of cores of their mesh
S=../common/submit_array.sh
for spec in "1 12:00:00 4G 72" "2 12:00:00 8G 36" "4 16:00:00 16G 36" "8 24:00:00 32G 36" \
            "16 2-00:00:00 64G 2" "32 2-00:00:00 96G 1"; do
    set -- $spec
    N=$(python -m farm dfg2d2_hpc count nprocs=$1 2>/dev/null || echo "$4")
    [ "$N" = "0" ] && continue
    J=$(bash $S dfg2d2_hpc "$HERE" "j2-np$1" "$1" "$2" "$3" "$N" "$SPIN" "nprocs=$1")
    echo "runs on $1 core(s): $N tasks, array job $J"
done
echo "monitor: squeue -u \$USER;  python -m farm dfg2d2_hpc status  (after: source ../config.sh; activate_python)"
