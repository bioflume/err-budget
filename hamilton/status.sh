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
