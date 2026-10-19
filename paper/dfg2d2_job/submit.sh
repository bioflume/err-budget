#!/bin/bash
# Submit the spin-up array, then the task array once every spin-up has finished.
set -e
cd "$(dirname "$0")"
source ./env.sh
mkdir -p logs
SPIN=$(sbatch --parsable slurm_spin.sh)
TASK=$(sbatch --parsable --dependency=afterok:$SPIN slurm_tasks.sh)
echo "spin-up array $SPIN, task array $TASK (starts when all spin-ups succeed)"
echo "monitor with: ./monitor.sh"
