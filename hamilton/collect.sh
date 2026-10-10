#!/bin/bash
# Run the analyses and pack the results to send back (a few MB).
cd "$(dirname "$0")"; HERE=$PWD
source ./config.sh; activate_python
mkdir -p results
(cd J1_dfg2d1 && python analyse_j1.py > ../results/j1.txt 2>&1)
(cd J2_dfg2d2 && python analyse_j2.py > ../results/j2.txt 2>&1)
(cd J3_qe && python analyse_j3.py > ../results/j3.txt 2>&1)
(cd J4_parallel && python analyse_j4.py > ../results/j4.txt 2>&1)
(cd J5_abc3d && python analyse_j5.py > ../results/j5.txt 2>&1)
# the raw per-task results (small), so that everything can be re-analysed
(cd "$WORK" && find . -path ./smoke -prune -o \( -name result.json -o -name series.npz \
    -o -name tensor.pkl -o -name ref.json -o -name error.txt \) -print \
    | tar czf "$HERE/results/raw.tar.gz" -T -)
tar czf hamilton_results.tar.gz results
ls -la hamilton_results.tar.gz
tail -n 3 results/j*.txt
