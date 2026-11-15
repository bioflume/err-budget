#!/bin/bash
# One task of a sweep per array element:  driver.py --sweep $SWEEP --index $SLURM_ARRAY_TASK_ID
# Submitted by slurm/submit.sh, which sets --array, --time and --mem per resource class.
# Ten minutes before the time limit Slurm sends USR1 to this script, which passes it to the
# driver: the driver stops after the run in progress, and the next submission goes on.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --signal=B:USR1@600
set -o pipefail
: "${SWEEP:?SWEEP not set}" "${REUSE_CODE:?REUSE_CODE not set}"
# shellcheck source=/dev/null
source "$REUSE_CODE/slurm/env.sh" || exit 1
cd "$REUSE_CODE" || exit 1
echo "$(date '+%Y-%m-%dT%H:%M:%S') $(hostname) sweep=$SWEEP index=$SLURM_ARRAY_TASK_ID"
python3 driver.py --sweep "$SWEEP" --index "$SLURM_ARRAY_TASK_ID" &
pid=$!
on_usr1() {
    echo "$(date '+%Y-%m-%dT%H:%M:%S') USR1: time limit in 10 minutes, stopping after the current run"
    kill -USR1 "$pid" 2>/dev/null
}
trap on_usr1 USR1
rc=0
while true; do
    wait "$pid"; rc=$?
    kill -0 "$pid" 2>/dev/null || break        # wait returned because of the trap: keep waiting
done
echo "$(date '+%Y-%m-%dT%H:%M:%S') driver exit $rc (0 done, 2 error recorded, 3 interrupted: resubmit, 4 gave up, 5 running in another job)"
exit "$rc"
