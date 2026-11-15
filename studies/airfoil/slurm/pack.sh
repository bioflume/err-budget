#!/bin/bash
# Pack a stage for the laptop (RUNBOOK.md):  bash slurm/pack.sh stage0|stage1|stage2|live
# Writes $AIRFOIL_WORK/packs/airfoil-<stage>-<date>.tar.gz with every stored run (JSON and
# compact log; the store is small and cumulative), the error files, the Slurm logs of
# the stage, sacct's records of its jobs, the sweep manifests, the configuration, the
# analysis output and the environment record; then prints the copy line for the Mac.
# Packed while jobs run, the archive is a snapshot (GNU tar may exit 1: "file changed").
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
stage=${1:?usage: pack.sh diag|numerics|convmap|confirm|settle|stage0|stage1|stage2|live}
# shellcheck source=/dev/null
source "$here/env.sh" python >/dev/null || exit 1
cd "$AIRFOIL_CODE" || exit 1
W=$AIRFOIL_WORK
stamp=$(date +%Y%m%d-%H%M)
meta="$W/packs/meta-$stage-$stamp"
mkdir -p "$meta"
case "$stage" in
    diag) sweeps="diag" ;;
    numerics) sweeps="numerics" ;;
    convmap) sweeps="convmap" ;;
    confirm) sweeps="confirm" ;;
    settle) sweeps="settle" ;;
    stage0) sweeps="check stage0" ;;
    stage1) sweeps="stage1 stage1ref" ;;
    stage2) sweeps="stage2" ;;
    live) sweeps="" ;;
    *) echo "unknown stage $stage" >&2; exit 1 ;;
esac
for s in $sweeps; do
    python3 driver.py --sweep "$s" --list --json > "$meta/manifest_$s.jsonl" 2>/dev/null
    python3 driver.py --sweep "$s" --status >> "$meta/status.txt" 2>&1
done
if [ -f "$W/jobs.txt" ]; then
    cp "$W/jobs.txt" "$meta/"
    ids=$(awk '{print $2}' "$W/jobs.txt" | paste -sd, -)
    sacct -j "$ids" -P --format=JobID,JobName%30,Partition,AllocCPUS,ReqMem,State,ExitCode,Elapsed,CPUTimeRAW,MaxRSS,NodeList,Start,End \
        > "$meta/sacct.txt" 2>&1
fi
cp -r "$AIRFOIL_CODE/config" "$meta/config"
for f in check_environment.txt module_avail.txt tutorial_airFoil2D_fvSolution; do
    [ -f "$W/$f" ] && cp "$W/$f" "$meta/"
done
(cd "$AIRFOIL_CODE" && git rev-parse HEAD 2>/dev/null || cat VERSION 2>/dev/null) > "$meta/code_version.txt"
out="$W/packs/airfoil-$stage-$stamp.tar.gz"
items=("runs" "packs/meta-$stage-$stamp")
[ -d "$W/results" ] && items+=("results")
[ -d "$W/live" ] && items+=("live")
logs=$(cd "$W" && ls logs/ 2>/dev/null | grep -E "^(af-($(echo "$sweeps" | tr ' ' '|'))_|check_|live_)" | sed 's|^|logs/|')
# shellcheck disable=SC2086
tar czf "$out" -C "$W" "${items[@]}" $logs
rc=$?
[ $rc -eq 1 ] && echo "note: tar exit 1 (a file changed while it was read): this archive is a snapshot"
[ $rc -gt 1 ] && { echo "tar failed ($rc)"; exit $rc; }
echo "packed $(du -h "$out" | cut -f1)  $out"
echo
echo "On the Mac (a Terminal there, not this session):"
echo "  scp grbj73@hamilton8.dur.ac.uk:$out \"\$AF/results/\""
