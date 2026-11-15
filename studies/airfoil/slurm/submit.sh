#!/bin/bash
# Submit a sweep (RUNBOOK.md):
#   bash slurm/submit.sh check                 # the 15-minute check, test partition
#   bash slurm/submit.sh stage0|stage1|stage1ref|stage2 [--dry-run]
#   bash slurm/submit.sh live                  # optional live errbudget run
#   AIRFOIL_MEM=8G AIRFOIL_TIME=1-00:00:00 bash slurm/submit.sh stage1   # other limits (OOM, time-out)
# For a sweep it submits only the runs that have no result yet, as one Slurm array per
# resource class (time limit, memory, cores; driver.py --groups), so the same command
# resubmits failed, interrupted or timed-out runs.  Every job id is appended to
# $AIRFOIL_WORK/jobs.txt for status.sh and pack.sh.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
sweep=${1:?usage: submit.sh check|diag|numerics|convmap|confirm|settle|stage0|stage1|stage1ref|stage2|live [--dry-run]}
dry=${2:-}
# shellcheck source=/dev/null
source "$here/env.sh" python || exit 1
set -eo pipefail
cd "$AIRFOIL_CODE"
log="$AIRFOIL_WORK/logs"
record() { echo "$(date '+%Y-%m-%dT%H:%M:%S') $1 $sweep $2" >> "$AIRFOIL_WORK/jobs.txt"; }

run() {
    echo "+ $*"
    if [ "$dry" = "--dry-run" ]; then return 0; fi
    jid=$("$@")
    echo "  submitted job $jid"
    record "$jid" "${desc:-}"
}

case "$sweep" in
    check)
        desc=check run sbatch --parsable $AIRFOIL_SB_EXTRA -p "$AIRFOIL_TEST_PARTITION" \
            -o "$log/check_%j.out" --export=ALL,AIRFOIL_CODE="$AIRFOIL_CODE" "$here/check.slurm"
        echo "when it has run:  tail -n 25 $log/check_<job>.out    (last line: CHECK PASSED)"
        ;;
    live)
        desc=live run sbatch --parsable $AIRFOIL_SB_EXTRA -p "$AIRFOIL_PARTITION" \
            -o "$log/live_%j.out" --export=ALL,AIRFOIL_CODE="$AIRFOIL_CODE" "$here/live.slurm"
        ;;
    diag|numerics|convmap|confirm|settle|stage0|stage1|stage1ref|stage2)
        groups=$(python3 driver.py --sweep "$sweep" --groups)
        if [ -z "$groups" ]; then
            echo "$sweep: every run has a result (or was given up); nothing to submit"
            python3 driver.py --sweep "$sweep" --status
            exit 0
        fi
        while read -r t mem np arr; do
            # after an OUT_OF_MEMORY or TIMEOUT: AIRFOIL_MEM=8G / AIRFOIL_TIME=1-00:00:00 bash slurm/submit.sh ...
            mem=${AIRFOIL_MEM:-$mem}; t=${AIRFOIL_TIME:-$t}
            desc="time=$t mem=$mem np=$np array=$arr"
            run sbatch --parsable $AIRFOIL_SB_EXTRA -p "$AIRFOIL_PARTITION" -J "af-$sweep" \
                --array="$arr" --time="$t" --mem="$mem" --cpus-per-task="$np" \
                -o "$log/%x_%A_%a.out" \
                --export=ALL,SWEEP="$sweep",AIRFOIL_CODE="$AIRFOIL_CODE" "$here/run_array.slurm"
        done <<< "$groups"
        python3 driver.py --sweep "$sweep" --status
        ;;
    *)
        echo "unknown sweep $sweep" >&2; exit 1 ;;
esac
