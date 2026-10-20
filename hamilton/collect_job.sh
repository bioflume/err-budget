#!/bin/bash
# Run collect.sh as a batch job: the J1 and J2 analyses replay the program over their
# grids, which takes 1-2 hours on one core (too long for a login node).
#   bash collect_job.sh        # then: ls -la hamilton_results.tar.gz when the job ends
set -e
cd "$(dirname "$0")"; HERE=$PWD
source ./config.sh
sbatch $SB_EXTRA -p "$SB_SHARED" -n 1 -t 08:00:00 --mem=16G -J errbudget-collect \
    -o "$HERE/collect_%j.out" --wrap "bash $HERE/collect.sh"
