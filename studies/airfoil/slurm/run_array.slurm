#!/bin/bash
# One run of a sweep per array task:  driver.py --sweep $SWEEP --index $SLURM_ARRAY_TASK_ID
# Submitted by slurm/submit.sh, which sets --array, --time, --mem, --cpus-per-task and
# the output file per resource class.  Ten minutes before the time limit Slurm sends
# USR1 to this script, which passes it to the driver: OpenFOAM writes its fields and
# stops, and the next submission resumes the run from there.
#SBATCH --ntasks=1
#SBATCH --signal=B:USR1@600
set -o pipefail
: "${SWEEP:?SWEEP not set}" "${AIRFOIL_CODE:?AIRFOIL_CODE not set}"
# shellcheck source=/dev/null
source "$AIRFOIL_CODE/slurm/env.sh" || exit 1
cd "$AIRFOIL_CODE" || exit 1
echo "$(date '+%Y-%m-%dT%H:%M:%S') $(hostname) sweep=$SWEEP index=$SLURM_ARRAY_TASK_ID cpus=${SLURM_CPUS_PER_TASK:-1}"
python3 driver.py --sweep "$SWEEP" --index "$SLURM_ARRAY_TASK_ID" &
pid=$!
on_usr1() {
    echo "$(date '+%Y-%m-%dT%H:%M:%S') USR1: time limit in 10 minutes, stopping the run cleanly"
    kill -USR1 "$pid" 2>/dev/null
}
trap on_usr1 USR1
rc=0
while true; do
    wait "$pid"; rc=$?
    kill -0 "$pid" 2>/dev/null || break        # wait returned because of the trap: keep waiting
done
echo "$(date '+%Y-%m-%dT%H:%M:%S') driver exit $rc (0 result stored, 2 error recorded, 3 interrupted: resubmit, 4 gave up, 5 running in another job)"
exit "$rc"
