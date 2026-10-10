#!/bin/bash
# One grid setting (or reference run) per array element; 243 tasks (0..242).
# Tasks are ordered cheapest first; the last three are the reference runs and
# take the longest (several hours each).  A finished task is skipped, so the
# array can simply be resubmitted after a failure or a time-out.
#SBATCH --job-name=dfg2-task
#SBATCH --array=0-242%40
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/task_%a.out
set -e
# SLURM runs a copy of this script, so use the submission directory
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source ./env.sh
cd "$KNOBS_ROOT/paper"
$PYTHON exp56_dfg2d2_build.py task $SLURM_ARRAY_TASK_ID
