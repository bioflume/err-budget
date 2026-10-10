#!/bin/bash
# submit_array.sh MODULE JOBDIR NAME NPROCS TIME MEM NARRAY [DEPENDENCY] [FILTER ...]
#
# Submits a SLURM array of NARRAY elements; element K runs the farm tasks
# K, K+NARRAY, ... of MODULE that match the FILTERs (key=value), each with
# NPROCS cores.  Prints the job id.  Used by the jobs' submit.sh scripts.
set -e
MOD=$1; JOBDIR=$2; NAME=$3; NP=$4; TIME=$5; MEM=$6; NA=$7; DEP=${8:-}
shift 8 || shift $#
FILTERS="$*"
source "$(dirname "$0")/../config.sh"
mkdir -p "$JOBDIR/logs"
D=""
[ -n "$DEP" ] && D="--dependency=afterok:$DEP"
sbatch --parsable $SB_EXTRA $D -p "$SB_SHARED" -J "$NAME" -N 1 -n "$NP" -c 1 --mem="$MEM" -t "$TIME" \
    --array=0-$((NA - 1)) -o "$JOBDIR/logs/${NAME}_%a.out" <<EOF
#!/bin/bash
source "$KNOBS_ROOT/hamilton/config.sh"
activate_python
load_openfoam
cd "$JOBDIR"
python -m farm $MOD array \$SLURM_ARRAY_TASK_ID $NA $FILTERS
EOF
