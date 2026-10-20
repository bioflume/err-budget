#!/bin/bash
# Progress of the build, the queue, and the last lines of any failed task.
cd "$(dirname "$0")"
source ./env.sh
cd "$ERRBUDGET_ROOT/paper"
$PYTHON exp56_dfg2d2_build.py status
if command -v squeue >/dev/null; then
  echo; squeue -u "$USER" -n dfg2-spin,dfg2-task -o "%.18i %.10j %.8T %.10M %R" | head -20
fi
echo; echo "tasks whose log reports an error:"
grep -l -i -E "error|traceback" "$ERRBUDGET_ROOT"/paper/dfg2d2_job/logs/*.out 2>/dev/null | head
