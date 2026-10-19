#!/bin/bash
# Run the whole build on one machine with N parallel workers (default: all cores).
# Safe to stop and restart at any time: finished tasks are kept, and an
# interrupted task is simply rerun.
#     ./run_local.sh 8          # 8 workers
#     ./monitor.sh              # progress
cd "$(dirname "$0")"
source ./env.sh
N=${1:-$(nproc)}
mkdir -p logs
cd "$KNOBS_ROOT/paper"
for w in $(seq 1 $N); do
  nohup $PYTHON exp56_dfg2d2_build.py worker >> dfg2d2_job/logs/worker_$w.out 2>&1 &
done
echo "$N workers started; logs in dfg2d2_job/logs/"
