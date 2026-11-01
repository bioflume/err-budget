#!/bin/bash
# Run the analyses and pack the results to send back (a few MB).
# Revision 4: the J1 and J2 analyses replay the program over their grids (1-2 hours on
# one core); run this as a batch job with  bash collect_job.sh  rather than on a login node.
# JOBS selects the analyses (default all), e.g.  JOBS="2 4" bash collect_job.sh  after the
# November 2026 reruns, whose J1, J3 and J5 results are already final.
cd "$(dirname "$0")"; HERE=$PWD
source ./config.sh; activate_python
mkdir -p results
JOBS=${JOBS:-"1 2 3 4 5"}
for j in $JOBS; do
    case $j in
        1) (cd J1_dfg2d1 && python analyse_j1.py > ../results/j1.txt 2>&1) ;;
        2) (cd J2_dfg2d2 && python analyse_j2.py > ../results/j2.txt 2>&1) ;;
        3) (cd J3_qe && python analyse_j3.py > ../results/j3.txt 2>&1) ;;
        4) (cd J4_parallel && python analyse_j4.py > ../results/j4.txt 2>&1) ;;
        5) (cd J5_abc3d && python analyse_j5.py > ../results/j5.txt 2>&1) ;;
    esac
done
# the raw per-task results (small), so that everything can be re-analysed; results that
# were moved aside (before_*, failed_*) and the smoke tests are left out
(cd "$WORK" && find . \( -path ./smoke -o -path './before_*' -o -path './failed_*' \) -prune -o \
    \( -name result.json -o -name series.npz \
    -o -name tensor.pkl -o -name ref.json -o -name error.txt -o -name residuals.npz \
    -o -name wall_geometry.json -o -name SPUN_UP -o -name SPIN_FAILED -o -name INLET_ANALYTIC \
    -o -name MAPPED -o -name attempts.txt \) -print \
    | tar czf "$HERE/results/raw.tar.gz" -T -)
tar czf hamilton_results.tar.gz results
ls -la hamilton_results.tar.gz
tail -n 3 results/j*.txt
