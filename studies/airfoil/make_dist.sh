#!/bin/bash
# Build dist/airfoil-batch.tar.gz, everything Hamilton needs (run on the Mac):
#   bash make_dist.sh
# The tracked files of this directory except tests/ and dist/, plus a copy of the
# errbudget engine (src/errbudget of the article's repository) as engine/errbudget for
# the live run, plus VERSION (commit, engine origin and hashes).  Writes dist/SHA256SUMS.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
engine_src=${ERRBUDGET_SRC:-$HOME/OpenResearch/knobs-paper/knobs-paper/src}
[ -f "$engine_src/errbudget/_core.py" ] || { echo "no errbudget engine at $engine_src" >&2; exit 1; }
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
cd "$here"
git ls-files . | grep -v -E '^(tests/|dist/)' | while read -r f; do
    mkdir -p "$stage/$(dirname "$f")"; cp -p "$f" "$stage/$f"
done
mkdir -p "$stage/engine"
rsync -a --exclude '__pycache__' --exclude '.DS_Store' "$engine_src/errbudget" "$stage/engine/"
commit=$(git rev-parse HEAD)
dirty=$( (git status --porcelain -- . | grep -v -E '^\?\? |dist/' || true) | wc -l | tr -d ' ')
engine_commit=$(git -C "$engine_src/.." rev-parse HEAD 2>/dev/null || echo unknown)
{
    echo "airfoil study, branch study/airfoil, commit $commit (uncommitted changes: $dirty)"
    echo "built $(date '+%Y-%m-%d %H:%M') on $(hostname)"
    echo "errbudget engine from $engine_src (repository commit $engine_commit):"
    (cd "$stage/engine" && find errbudget -name '*.py' | sort | xargs shasum -a 256)
} > "$stage/VERSION"
mkdir -p dist
# no macOS extended attributes: GNU tar on Hamilton would warn about each file
COPYFILE_DISABLE=1 tar --no-xattrs --no-mac-metadata -czf dist/airfoil-batch.tar.gz -C "$stage" .
(cd dist && shasum -a 256 airfoil-batch.tar.gz > SHA256SUMS)
echo "built dist/airfoil-batch.tar.gz ($(du -h dist/airfoil-batch.tar.gz | cut -f1)), $(tar tzf dist/airfoil-batch.tar.gz | wc -l | tr -d ' ') entries"
cat dist/SHA256SUMS
