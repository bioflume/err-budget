#!/bin/bash
# Optional (RUNBOOK.md, step 10): errbudget itself, live, through its job file and the
# command interface, from the user's default start at the hardest design.  After
# Stage 1 every run it asks for is already stored, so it takes minutes and must return
# the setting the replay returns.  Submitted by:  bash slurm/submit.sh live
#SBATCH --job-name=af-live
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=12:00:00
#SBATCH --signal=B:USR1@600
set -o pipefail
# shellcheck source=/dev/null
source "${AIRFOIL_CODE:?AIRFOIL_CODE not set}/slurm/env.sh" || exit 1
cd "$AIRFOIL_CODE" || exit 1
mkdir -p "$AIRFOIL_WORK/live"
python3 live_tune.py errbudget_job.yaml --out "$AIRFOIL_WORK/live"
