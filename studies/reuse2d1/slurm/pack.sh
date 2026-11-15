#!/bin/bash
# Pack the results for the Mac (RUNBOOK.md, step 6):  bash slurm/pack.sh
# Writes $REUSE_WORK/packs/reuse2d1-<date>.tar.gz with every stored run, the task error
# files, the Slurm logs, sacct's records, the task lists and the environment record, and
# prints the copy line for the Mac.  Packed while jobs run, the archive is a snapshot.
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=/dev/null
source "$here/env.sh" python >/dev/null || exit 1
cd "$REUSE_CODE" || exit 1
W=$REUSE_WORK
stamp=$(date +%Y%m%d-%H%M)
meta="$W/packs/meta-$stamp"
mkdir -p "$meta"
for s in check grid ref; do
    python3 driver.py --sweep "$s" --list --json > "$meta/tasks_$s.jsonl" 2>/dev/null
    python3 driver.py --sweep "$s" --status >> "$meta/status.txt" 2>&1
done
if [ -f "$W/jobs.txt" ]; then
    cp "$W/jobs.txt" "$meta/"
    ids=$(awk '{print $2}' "$W/jobs.txt" | paste -sd, -)
    sacct -j "$ids" -P --format=JobID,JobName%30,Partition,AllocCPUS,ReqMem,State,ExitCode,Elapsed,CPUTimeRAW,MaxRSS,NodeList,Start,End \
        > "$meta/sacct.txt" 2>&1
fi
[ -f "$W/check_environment.txt" ] && cp "$W/check_environment.txt" "$meta/"
cat "$REUSE_CODE/VERSION" > "$meta/code_version.txt" 2>/dev/null
out="$W/packs/reuse2d1-$stamp.tar.gz"
tar czf "$out" -C "$W" runs tasks logs "packs/meta-$stamp"
rc=$?
[ $rc -eq 1 ] && echo "note: tar exit 1 (a file changed while it was read): this archive is a snapshot"
[ $rc -gt 1 ] && { echo "tar failed ($rc)"; exit $rc; }
echo "packed $(du -h "$out" | cut -f1)  $out"
echo
echo "On the Mac (a Terminal there, not this session):"
echo "  scp grbj73@hamilton8.dur.ac.uk:$out \"\$RE/results/\""
