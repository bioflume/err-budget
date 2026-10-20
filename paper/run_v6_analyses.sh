#!/usr/bin/env bash
# Revision 4 analyses on the article's machine (the one with ABINIT 9.10.4, LAMMPS
# 20240207 and the GROMACS builds that produced the tensors in data/).
#
#   bash run_v6_analyses.sh check                      # environment only (a minute)
#   bash run_v6_analyses.sh before                     # everything that does not need Hamilton
#   bash run_v6_analyses.sh after hamilton_results.tar.gz   # once the cluster results are back
#   bash run_v6_analyses.sh collect                    # RESULTS_SUMMARY.md and the archive
#
# Every step writes its output to ../logs/v6/<step>.txt and its exit status to
# ../logs/v6/STATUS.txt; a failed step is reported and the script goes on, except the
# environment check and the tests.  KNOBS_JOBS (default 2) pooled configurations run at once.
set -u -o pipefail
ORIG=$PWD
cd "$(dirname "$0")"
PAPER=$PWD
LOG=$PAPER/../logs/v6
mkdir -p "$LOG" data figures
JOBS=${KNOBS_JOBS:-2}
PY=${PYTHON:-python}
export PYTHONPATH="$PAPER/../src:$PAPER:${PYTHONPATH:-}"
FAILED=()

step() {        # step NAME command ...
    local name=$1; shift
    echo "== $(date '+%H:%M:%S') $name: $*"
    "$@" > "$LOG/$name.txt" 2>&1
    local st=$?
    echo "$(date '+%F %T') $name exit $st" >> "$LOG/STATUS.txt"
    if [ $st -ne 0 ]; then
        echo "   FAILED ($st): see logs/v6/$name.txt"; tail -5 "$LOG/$name.txt"; FAILED+=("$name")
    else
        tail -2 "$LOG/$name.txt"
    fi
    return $st
}

pooled_parallel() {   # pooled_parallel TAG [--only KEYS] CONFIG ...
    local tag=$1; shift
    local only=()
    if [ "${1:-}" = "--only" ]; then only=(--only "$2"); shift 2; fi
    printf '%s\n' "$@" | xargs -P "$JOBS" -I{} bash -c \
        "$PY pooled.py ${only[*]:-} {} > '$LOG/${tag}_{}.txt' 2>&1 && echo '   done {}' || echo '   FAILED {}'"
    local bad
    bad=$(grep -L "^TOTAL" "$LOG/${tag}_"*.txt 2>/dev/null || true)
    [ -z "$bad" ] || { echo "   pooled runs without a TOTAL line: $bad"; FAILED+=("$tag"); }
    echo "$(date '+%F %T') $tag pooled runs: $*" >> "$LOG/STATUS.txt"
}

MAIN="v1.2.1 v1.2 v1.1 v1.1+gated"
ABL121="v1.2.1-nojoint v1.2.1-nolocal v1.2.1-credit v1.2.1-constant v1.2.1-nores \
        v1.2.1-accept-unconverged v1.2.1-nested-as-disc v1.2.1-screened"
ABL12="v1.2-nojoint v1.2-nolocal v1.2-credit v1.2-constant v1.2-nores \
       v1.2-accept-unconverged v1.2-nested-as-disc v1.2-screened"

check() {
    $PY - <<'EOF' || exit 1
import importlib, shutil, sys
ok = True
for m in ("numpy", "scipy", "matplotlib", "sympy", "pytest", "knobs"):
    try:
        mod = importlib.import_module(m)
        print(f"{m:<10} {getattr(mod, '__version__', '')} {getattr(mod, '__file__', '')}")
    except Exception as exc:
        print(f"{m:<10} MISSING ({exc})"); ok = False
try:
    import lammps; print("lammps     python module present")
except Exception as exc:
    print(f"lammps     python module MISSING ({exc}): exp74 cannot run")
for exe in ("abinit", "gmx", "gmx_d", "git"):
    print(f"{exe:<10} {shutil.which(exe) or 'MISSING'}")
import knobs
print("knobs version", knobs.__version__, "V121" if hasattr(knobs, "V121") else "NO V121")
sys.exit(0 if ok and hasattr(knobs, "V121") else 1)
EOF
    git -C .. tag | grep -q v1.2.0-frozen && echo "git tag v1.2.0-frozen present" \
        || { [ -d frozen_src/knobs ] && echo "no git tag; paper/frozen_src present" \
             || { echo "neither the git tag v1.2.0-frozen nor paper/frozen_src: held-out replays cannot run"; exit 1; }; }
}

before() {
    check || { echo "environment check failed"; exit 1; }
    step A01_tests $PY -m pytest ../tests -q || { echo "tests failed: stop"; exit 1; }
    step A02_spec_tests $PY ../tests/test_spec_v121.py || { echo "specification tests failed: stop"; exit 1; }
    # the ABINIT reference first: its uncertainty decides which ABINIT targets are graded
    step A03_abinit_reference $PY exp75_abinit_reference.py
    grep -q "REFERENCE REPRODUCED" "$LOG/A03_abinit_reference.txt" \
        || echo "   WARNING: the ABINIT reference was not reproduced; cases.py keeps the revision-2 uncertainty"
    step A04_pppm_mechanism $PY exp74_pppm_mechanism.py
    echo "== $(date '+%H:%M:%S') pooled tables: $MAIN $ABL121 $ABL12 ($JOBS at a time)"
    pooled_parallel A05_pooled $MAIN $ABL121 $ABL12
    step A06_summary $PY summ.py v1.2.1 v1.2 v1.1
    step A07_baselines env KNOBS_CFG=v1.2.1 $PY exp72_baselines_pooled.py
    step A08_baselines_lmp env KNOBS_CFG=v1.2.1 KNOBS_BASE_OUT=baselines_lmp_gmx.json \
        $PY exp72_baselines_pooled.py LMP-20 LMP-40 LMP-64 GMX-d
    step A09_accounting $PY exp76_accounting.py
    step A10_reltol_2d1 $PY exp77_reltol_2d1.py
    step A11_gromacs_floor $PY exp78_gromacs_floor_v121.py
    step A12_floor_sensitivity $PY exp58_floor_sensitivity.py v1.2.1
    step A13_unconverged $PY exp79_unconverged.py
    step A14_decline_reasons $PY exp73_decline_reasons.py
    step A15_heldout_frozen $PY -W ignore exp62_heldout.py
    step A16_heldout_current $PY -W ignore exp62_heldout.py --current
    step A17_dfg2_detail $PY exp71_dfg2_detail.py
    step A18_misdeclare $PY exp65_misdeclare_v12.py
    step A19_estimate $PY exp66_estimate_v12.py
    step A20_margin $PY exp67_margin_v12.py
    step A21_transfer_floor $PY exp64_transfer_floor.py
    step A21b_numbers $PY exp80_v6_numbers.py
    step A21c_bound_check $PY exp81_bound_check.py
    step A22_figures $PY exp63_figures.py
    step A23_tables $PY exp70_tables.py
    report
}

after() {
    local arch=${1:?usage: run_v6_analyses.sh after hamilton_results.tar.gz}
    case "$arch" in /*) ;; *) arch="$ORIG/$arch" ;; esac
    [ -f "$arch" ] || { echo "no such archive: $arch"; exit 1; }
    step B01_hamilton $PY exp60_hamilton.py "$arch"
    if [ -f ../hamilton/results/j2_tensor.pkl ] && $PY -c "import cases; cases._dfg2_refined()" 2>/dev/null; then
        echo "== $(date '+%H:%M:%S') pooled tables, refined 2D-2 study only (merged into the files of 'before')"
        pooled_parallel B02_pooled_2d2h --only 2D-2h $MAIN $ABL121 $ABL12
        step B03_baselines_2d2h env KNOBS_CFG=v1.2.1 $PY exp72_baselines_pooled.py 2D-2h
    else
        echo "   no refined 2D-2 case (hamilton/results/j2_tensor.pkl with a case specification): B02-B03 skipped"
        FAILED+=("B02_no_2d2h")
    fi
    step B04_accounting $PY exp76_accounting.py
    step B05_unconverged $PY exp79_unconverged.py
    step B06_decline_reasons $PY exp73_decline_reasons.py
    # (B01 has replayed H1-H3 with the frozen program and with v1.2.1: exp60 runs exp62 twice)
    step B07_tables $PY exp70_tables.py
    report
}

collect() {
    step C01_collect $PY exp82_collect_v6.py
}

report() {
    if [ ${#FAILED[@]} -eq 0 ]; then echo "ALL STEPS FINISHED"; return 0; fi
    echo "FAILED STEPS: ${FAILED[*]}"
    return 1
}

case "${1:-}" in
    check) check ;;
    before) before ;;
    after) shift; after "$@" ;;
    collect) collect ;;
    *) sed -n 2,12p "$0"; exit 2 ;;
esac
