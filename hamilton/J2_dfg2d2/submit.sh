#!/bin/bash
# J2: DFG 2D-2 refined (about 1,500-2,000 core-hours with the spin-ups; 1-3 days).
#   bash submit.sh          # one spin-up job per mesh, coarse to fine, then the runs of each mesh
#   bash submit.sh smoke    # a 15-minute check (two tiny meshes, the finer one mapped and parallel)
#
# November 2026: each mesh is spun up by its own job, on the cores its runs use and with a
# one-day limit, after the job of the next coarser mesh has succeeded; the runs of a mesh
# start when its spin-up has succeeded, and SLURM cancels the jobs that depend on a failed
# one (--kill-on-invalid-dep).  Results and meshes made before the inlet fix are refused.
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ../config.sh
activate_python
mkdir -p logs
if [ "${1:-}" = "smoke" ]; then
    sbatch --parsable $SB_EXTRA -p test -N 1 -n 2 -t 00:15:00 --mem=4G -J j2-smoke \
        -o logs/smoke_%j.out --export=ALL,ERRBUDGET_SMOKE=1 <<EOF
#!/bin/bash
source "$ERRBUDGET_ROOT/hamilton/config.sh"; activate_python; load_openfoam; cd "$HERE"
# a 10-second window, so that the watchdog reads the logs of the container runs here too
export ERRBUDGET_SMOKE=1 WORK="$WORK/smoke" ERRBUDGET_OF_WATCH_WINDOW=10
rm -rf "\$WORK/dfg2d2_hpc" "\$WORK/dfg2d2_cases"   # an earlier smoke test would be skipped, not re-run
python dfg2d2_hpc.py spin && python -m farm dfg2d2_hpc array 0 1 && python -m farm dfg2d2_hpc status
cat \$WORK/dfg2d2_hpc/*/error.txt 2>/dev/null | tail -n 20
NOFAIL=\$(python -c "import glob,json; f=glob.glob('\$WORK/dfg2d2_hpc/*/result.json'); print(int(len(f) == 3 and all(json.load(open(x))['failed'] is None for x in f)))")
if [ -f "\$WORK/dfg2d2_cases/mesh_16/INLET_ANALYTIC" ] && [ "\$(ls \$WORK/dfg2d2_hpc/*/series.npz 2>/dev/null | wc -l)" = 3 ] \
   && [ "\$NOFAIL" = 1 ] && ! ls \$WORK/dfg2d2_hpc/*/error.txt > /dev/null 2>&1; then
    echo "SMOKE TEST PASSED"
else
    echo "SMOKE TEST FAILED"
fi
EOF
    exit 0
fi

# 0. refuse to submit twice, or over results and meshes from before the inlet fix
if squeue -h -u "$USER" -o "%j" 2>/dev/null | grep -q '^j2-'; then
    echo "J2 jobs are already queued or running (squeue -u $USER); submitting again would run tasks twice."
    exit 1
fi
python - <<'PYEOF' || exit 1
import glob, json, os, sys
W = os.environ["WORK"]
C = os.path.join(W, "dfg2d2_cases")
bad = [os.path.dirname(f) for f in glob.glob(os.path.join(C, "mesh_*", "MAPPED"))
       if not os.path.exists(os.path.join(os.path.dirname(f), "INLET_ANALYTIC"))]
old = [f for f in glob.glob(os.path.join(W, "dfg2d2_hpc", "*", "result.json"))
       if json.load(open(f)).get("inlet") != "analytic"]
if bad or old:
    msg = []
    if bad:
        msg.append(f"{len(bad)} meshes were mapped before the inlet fix: {' '.join(sorted(bad))}")
    if old:
        msg.append(f"{len(old)} results in {W}/dfg2d2_hpc come from before the inlet fix")
    sys.exit("\n".join(msg) + "\nMove them aside (do not delete them), e.g.\n"
             f"  mkdir -p {W}/before_inlet_fix\n"
             f"  mv {W}/dfg2d2_hpc {W}/before_inlet_fix/\n"
             + "".join(f"  mv {d} {W}/before_inlet_fix/\n" for d in sorted(bad))
             + "and submit again.")
PYEOF

[ -f "$WORK/J2_SUBMITTED" ] || date +%F > "$WORK/J2_SUBMITTED"   # for the core-hour guard (status.sh)

# 1. spin-ups: one job per mesh that is not spun up yet, each after the coarser one
declare -A SPINJOB
PREV=""
for spec in "24 1 4G" "32 1 4G" "48 2 8G" "64 4 16G" "96 8 32G" "128 16 64G" "192 32 96G"; do
    set -- $spec
    NA=$1; NP=$2; MEM=$3
    if python -c "import sys, dfg2d2_hpc as J; sys.exit(0 if J.spun($NA) else 1)" 2>/dev/null; then
        echo "n_arc $NA: already spun up (analytic inlet)"
        continue
    fi
    D=""
    [ -n "$PREV" ] && D="--dependency=afterok:$PREV --kill-on-invalid-dep=yes"
    J=$(sbatch --parsable $SB_EXTRA $D -p "$SB_SHARED" -N 1 -n "$NP" -c 1 -t 1-00:00:00 --mem="$MEM" \
        -J "j2-spin$NA" -o "logs/spin${NA}_%j.out" <<EOF
#!/bin/bash
source "$ERRBUDGET_ROOT/hamilton/config.sh"; activate_python; load_openfoam; cd "$HERE"
python dfg2d2_hpc.py spin $NA
EOF
)
    SPINJOB[$NA]=$J
    PREV=$J
    echo "n_arc $NA: spin-up job $J ($NP cores)"
done

# 2. the runs of each mesh, after its spin-up (references at n_arc 128 and 192 included)
S=../common/submit_array.sh
for spec in "24 1 12:00:00 4G" "32 1 12:00:00 4G" "48 2 12:00:00 8G" "64 4 16:00:00 16G" \
            "96 8 1-00:00:00 32G" "128 16 2-00:00:00 64G" "192 32 2-00:00:00 96G"; do
    set -- $spec
    NA=$1; NP=$2; TIME=$3; MEM=$4
    N=$(python -m farm dfg2d2_hpc count n_arc=$NA)
    if [ "$(python -m farm dfg2d2_hpc left n_arc=$NA)" = "0" ]; then
        echo "n_arc $NA: all $N tasks done"
        continue
    fi
    J=$(bash $S dfg2d2_hpc "$HERE" "j2-na$NA" "$NP" "$TIME" "$MEM" "$N" "${SPINJOB[$NA]:-}" "n_arc=$NA")
    echo "n_arc $NA: $N tasks, array job $J${SPINJOB[$NA]:+ (after spin-up job ${SPINJOB[$NA]})}"
done
echo "monitor: squeue -u \$USER;  bash ../status.sh  (core-hours so far: bash ../coreh.sh j2)"
