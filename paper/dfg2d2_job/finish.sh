#!/bin/bash
# After every task has finished: assemble the grid and print the analysis.
cd "$(dirname "$0")"
source ./env.sh
cd "$KNOBS_ROOT/paper"
$PYTHON exp56_dfg2d2_build.py collect
$PYTHON exp56_dfg2d2.py | tee figures/sec8x_dfg2d2.txt
