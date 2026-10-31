#!/bin/bash
# Progress of every job: tasks done / with errors.
cd "$(dirname "$0")"
source ./config.sh; activate_python
(cd J1_dfg2d1 && python -m farm dfg2d1_hpc status)
(cd J2_dfg2d2 && python -m farm dfg2d2_hpc status)
(cd J3_qe && python -m farm qe_hpc status)
(cd J4_parallel && python -m farm parallel_hpc status)
(cd J5_abc3d && python -m farm abc_hpc status)
grep -l . $WORK/*/*/error.txt 2>/dev/null | head -20 | sed 's/^/error: /'
# core-hours of the November 2026 reruns, from their first submission (guard: J2 2,200)
[ -f "$WORK/J2_SUBMITTED" ] && { echo "J2 core-hours:"; bash coreh.sh j2 "$(cat "$WORK/J2_SUBMITTED")" 2200; }
[ -f "$WORK/J4_SUBMITTED" ] && { echo "J4 core-hours:"; bash coreh.sh j4 "$(cat "$WORK/J4_SUBMITTED")"; }
