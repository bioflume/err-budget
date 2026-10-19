#!/bin/bash
# Submit all five heavy-tier jobs.  See README.md.
set -e
cd "$(dirname "$0")"
for J in J1_dfg2d1 J2_dfg2d2 J3_qe J4_parallel J5_abc3d; do
    echo "== $J"; (cd $J && bash submit.sh "$@")
done
