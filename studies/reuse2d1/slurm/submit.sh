#!/bin/bash
# Submit (RUNBOOK.md):
#   bash slurm/submit.sh check             # the 15-minute check, test partition
#   bash slurm/submit.sh grid|ref [--dry-run]
#   REUSE_MEM=8G REUSE_TIME=1-00:00:00 bash slurm/submit.sh grid   # other limits (OOM, time-out)
# For a sweep it submits only the tasks that are not finished, not running and not given up,
# as one Slurm array per resource class, so the same command resubmits failed, interrupted
# or timed-out tasks.  Every job id is appended to $REUSE_WORK/jobs.txt.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
sweep=${1:?usage: submit.sh check|grid|ref [--dry-run]}
dry=${2:-}
# shellcheck source=/dev/null
source "$here/env.sh" python || exit 1
set -eo pipefail
cd "$REUSE_CODE"
log="$REUSE_WORK/logs"
record() { echo "$(date '+%Y-%m-%dT%H:%M:%S') $1 $sweep $2" >> "$REUSE_WORK/jobs.txt"; }
run() {
    echo "+ $*"
    if [ "$dry" = "--dry-run" ]; then return 0; fi
    jid=$("$@")
    echo "  submitted job $jid"
    record "$jid" "${desc:-}"
}
case "$sweep" in
    check)
        desc=check run sbatch --parsable $REUSE_SB_EXTRA -p "$REUSE_TEST_PARTITION" \
            -o "$log/check_%j.out" --export=ALL,REUSE_CODE="$REUSE_CODE" "$here/check.slurm"
        echo "when it has run:  tail -n 20 $log/check_<job>.out    (last line: CHECK PASSED)"
        ;;
    grid|ref)
        groups=$(python3 driver.py --sweep "$sweep" --groups)
        if [ -z "$groups" ]; then
            echo "$sweep: every task has its results, runs elsewhere or was given up; nothing to submit"
            python3 driver.py --sweep "$sweep" --status
            exit 0
        fi
        while read -r t mem arr; do
            mem=${REUSE_MEM:-$mem}; t=${REUSE_TIME:-$t}
            desc="time=$t mem=$mem array=$arr"
            run sbatch --parsable $REUSE_SB_EXTRA -p "$REUSE_PARTITION" -J "re-$sweep" \
                --array="$arr" --time="$t" --mem="$mem" \
                -o "$log/%x_%A_%a.out" \
                --export=ALL,SWEEP="$sweep",REUSE_CODE="$REUSE_CODE" "$here/run_array.slurm"
        done <<< "$groups"
        python3 driver.py --sweep "$sweep" --status
        ;;
    *)
        echo "unknown sweep $sweep" >&2; exit 1 ;;
esac
