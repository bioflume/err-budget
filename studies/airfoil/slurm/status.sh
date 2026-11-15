#!/bin/bash
# Progress and cost (RUNBOOK.md, step 7):  bash slurm/status.sh [sweep ...]
# Per sweep: runs done / failed (a result) / error / started / missing; then the jobs in
# the queue; then the core-hours used per stage from sacct, against the guards of
# PROTOCOL.md 12; then the newest error files.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$here/env.sh" python >/dev/null || exit 1
cd "$AIRFOIL_CODE" || exit 1
sweeps=${*:-check stage0 stage1 stage1ref stage2}
echo "== runs"
for s in $sweeps; do python3 driver.py --sweep "$s" --status 2>&1 | sed 's/^/  /'; done
echo "== queue"
squeue -u "$USER" -o "%.18i %.12j %.8T %.10M %.10l %.5C %R" 2>/dev/null | sed 's/^/  /'
echo "== core-hours (sacct; running jobs counted so far)"
if [ -f "$AIRFOIL_WORK/jobs.txt" ]; then
    for stage in check diag numerics convmap confirm settle stage0 stage1 stage2; do
        ids=$(awk -v s="$stage" '($3==s || (s=="stage1" && $3=="stage1ref")) {print $2}' "$AIRFOIL_WORK/jobs.txt" | paste -sd, -)
        [ -z "$ids" ] && continue
        sec=$(sacct -X -n -P -j "$ids" --format=CPUTimeRAW 2>/dev/null | awk '{s+=$1} END {print s+0}')
        guard_var="AIRFOIL_GUARD_$stage"
        guard=${!guard_var:-0}
        awk -v s="$sec" -v g="$guard" -v n="$stage" 'BEGIN {
            h = s / 3600; printf "  %-8s %8.1f core-hours", n, h;
            if (g > 0) { printf "  (guard %d)", g; if (h > g) printf "  GUARD PASSED: stop submitting %s and report", n }
            printf "\n" }'
    done
    echo "  states: $(sacct -X -n -P -j "$(awk '{print $2}' "$AIRFOIL_WORK/jobs.txt" | paste -sd, -)" --format=State 2>/dev/null | sort | uniq -c | tr -s ' ' | paste -sd, -)"
fi
echo "== newest error files"
ls -t "$AIRFOIL_WORK"/runs/*.error-*.json 2>/dev/null | head -5 | sed 's/^/  /'
