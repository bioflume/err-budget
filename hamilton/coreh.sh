#!/bin/bash
# Core-hours used by the jobs whose names start with PREFIX- (from sacct), since START.
#   bash coreh.sh j2 [START (YYYY-MM-DD, default 14 days ago)] [GUARD (core-hours)]
# With a GUARD it says so when the total passes it: stop submitting and report.
P=${1:?usage: coreh.sh PREFIX [START] [GUARD]}
START=${2:-$(date -d '14 days ago' +%F)}
GUARD=${3:-}
sacct -X -n -P -u "$USER" -S "$START" --format=JobName,CPUTimeRAW 2>/dev/null | \
awk -F'|' -v p="^$P-" -v lab="$P-*" -v g="$GUARD" -v s0="$START" '
$1 ~ p { s[$1] += $2; t += $2 }
END { for (k in s) printf "  %-22s %9.1f core-hours\n", k, s[k] / 3600
      printf "  %-22s %9.1f core-hours since %s\n", "total " lab, t / 3600, s0
      if (g != "" && t / 3600 > g) printf "GUARD PASSED: %.0f > %s core-hours: stop submitting and report\n", t / 3600, g }'
