#!/usr/bin/env bash
# Replays the program over the stored runs and regenerates every number, table and figure
# of the article and its supplement that needs no solver; README.md maps each result to
# its step.  The steps that run a solver are in reproduce.sh.
#
#   bash replay_all.sh               # every step below
#   bash replay_all.sh --no-eight    # without the eight-parameter study (6.5 CPU-hours)
#   bash replay_all.sh --pre         # also the PRE comparison, which downloads the PRE code
#
# Results go to data/ and figures/, where `git diff` compares them with the stored ones.
# The printed output of step NAME goes to ../logs/replay/NAME.txt.  ERRBUDGET_JOBS
# (default 2) pooled configurations run at once.
set -u
cd "$(dirname "$0")" || exit 1
LOG=../logs/replay
mkdir -p "$LOG"
PY=${PYTHON:-python}
JOBS=${ERRBUDGET_JOBS:-2}
EIGHT=1 PRE=0
for a in "$@"; do
    case $a in
        --no-eight) EIGHT=0 ;;
        --pre) PRE=1 ;;
        *) sed -n 2,12p "$0"; exit 2 ;;
    esac
done
FAILED=()

step() {
    local name=$1; shift
    echo "$name"
    "$@" > "$LOG/$name.txt" 2>&1 || { echo "  failed: see logs/replay/$name.txt"; FAILED+=("$name"); }
}

pooled() {
    printf '%s\n' "$@" | xargs -P "$JOBS" -I{} sh -c "$PY pooled.py {} > $LOG/pooled_{}.txt 2>&1" \
        || FAILED+=(pooled)
}

# The released program and its ablation (Tables 4, S-3.1, S-3.2)
step tests $PY -m pytest ../tests -q
echo "pooled (9 configurations, $JOBS at a time)"
pooled v1.2.1 v1.2.1-nojoint v1.2.1-nolocal v1.2.1-credit v1.2.1-constant v1.2.1-nores \
       v1.2.1-accept-unconverged v1.2.1-nested-as-disc v1.2.1-screened
step baselines env ERRBUDGET_CFG=v1.2.1 $PY exp72_baselines_pooled.py
step baselines_md env ERRBUDGET_CFG=v1.2.1 ERRBUDGET_BASE_OUT=baselines_lmp_gmx.json \
    $PY exp72_baselines_pooled.py LMP-20 LMP-40 LMP-64 GMX-d
step accounting $PY exp76_accounting.py
step decline_reasons $PY exp73_decline_reasons.py
step unconverged $PY exp79_unconverged.py
step estimate $PY exp66_estimate_v12.py
step separability $PY exp28_separability.py
step margin $PY exp67_margin_v12.py
step misdeclare $PY exp65_misdeclare_v12.py
step classify $PY exp17_classify.py
step transfer $PY exp64_transfer_floor.py
step section5 $PY exp40_section5.py
step bound_check $PY exp81_bound_check.py
step numbers $PY exp80_v6_numbers.py

# Pre-registered problems (Table 5, S-3.7)
step heldout_frozen $PY -W ignore exp62_heldout.py
step heldout_released $PY -W ignore exp62_heldout.py --current

# Production codes (Section 7, S-5)
step reltol_2d1 $PY exp77_reltol_2d1.py
step gromacs_floor $PY exp78_gromacs_floor_v121.py
step floor_sensitivity $PY exp58_floor_sensitivity.py v1.2.1
step dfg2d2_coarse $PY exp56_dfg2d2.py
step dfg2_detail $PY exp71_dfg2_detail.py
step dfg2_traces $PY -W ignore exp83_dfg2_traces.py
step dfg2_strict $PY -W ignore exp84_dfg2_strict.py

# Reruns with the released program (S-3.5, S-4.5)
step reuse_released $PY exp30_transfer.py
step reuse_hardest_released $PY exp37_calhardest.py
step spacing_released $PY exp34_spacing.py
step spacing_margin $PY exp34_spacing_margin.py
step spacing_estimate_released $PY exp35_whyfine.py
[ "$EIGHT" = 1 ] && step eight_released $PY exp54_eight.py

# Tables and figures
step tables $PY exp70_tables.py
step fig_2d1_l5_data $PY extract_2d1_l5.py
step fig_2d1_l5 $PY make_2d1_l5_figures.py
step fig_baselines $PY make_baselines_figure.py
step fig_lammps $PY make_lammps_figure.py
step fig_problems $PY make_solution_figure.py
step fig_pre $PY make_pre_figure.py

# Cluster studies, from their raw packs when these are unpacked (../DATA.md)
R=../studies/reuse2d1/results/reuse2d1-20261115-1111
if [ -d "$R/runs" ]; then
    step reuse2d1 $PY ../studies/reuse2d1/analysis/reuse.py --runs "$R/runs" --out "$R/res"
else
    echo "reuse2d1: raw pack not unpacked; stored summary in studies/reuse2d1/summary/"
fi
A=../studies/airfoil/results/stage1
if [ -d "$A/runs" ]; then
    step airfoil_stage1 $PY ../studies/airfoil/analysis/stage1.py --runs "$A/runs" \
        --config "$(ls -d "$A"/packs/meta-stage1-*/config)" --out "$A/analysis"
else
    echo "airfoil: raw pack not unpacked; stored summary in studies/airfoil/summary/"
fi

if [ "$PRE" = 1 ]; then
    P=../studies/pre_comparison
    step pre_fetch sh "$P/fetch_pre.sh"
    for s in arm1 check_replay arm2 summarize; do
        step "pre_$s" sh -c "cd $P && .venv/bin/python $s.py"
    done
fi

if [ ${#FAILED[@]} -eq 0 ]; then echo "all steps finished"; else echo "failed: ${FAILED[*]}"; exit 1; fi
