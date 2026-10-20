#!/bin/bash
# Spin up one mesh level per array element (indices 0..5 -> n_arc 12 16 24 32 48 64).
#SBATCH --job-name=dfg2-spin
#SBATCH --array=0-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=logs/spin_%a.out
set -e
# SLURM runs a copy of this script, so use the submission directory
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
source ./env.sh
LEVELS=(12 16 24 32 48 64)
cd "$ERRBUDGET_ROOT/paper"
$PYTHON exp56_dfg2d2_build.py spin ${LEVELS[$SLURM_ARRAY_TASK_ID]}
