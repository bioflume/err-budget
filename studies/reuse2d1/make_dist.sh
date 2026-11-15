#!/bin/bash
# Build dist/reuse2d1-batch.tar.gz, everything Hamilton needs (run on the Mac):
#   bash make_dist.sh
# The tracked files of this directory except tests/ and dist/, plus VERSION (commit).
# Writes dist/SHA256SUMS.  The errbudget engine is not needed on Hamilton: the search is
# replayed on the Mac from the stored runs (PROTOCOL.md 8).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
cd "$here"
git ls-files . | grep -v -E '^(tests/|dist/)' | while read -r f; do
    mkdir -p "$stage/$(dirname "$f")"; cp -p "$f" "$stage/$f"
done
commit=$(git rev-parse HEAD)
dirty=$( (git status --porcelain -- . | grep -v -E '^\?\? |dist/' || true) | wc -l | tr -d ' ')
{
    echo "reuse2d1 study, branch $(git branch --show-current), commit $commit (uncommitted changes: $dirty)"
    echo "built $(date '+%Y-%m-%d %H:%M') on $(hostname)"
} > "$stage/VERSION"
mkdir -p dist
COPYFILE_DISABLE=1 tar --no-xattrs --no-mac-metadata -czf dist/reuse2d1-batch.tar.gz -C "$stage" .
(cd dist && shasum -a 256 reuse2d1-batch.tar.gz > SHA256SUMS)
echo "built dist/reuse2d1-batch.tar.gz ($(du -h dist/reuse2d1-batch.tar.gz | cut -f1)), $(tar tzf dist/reuse2d1-batch.tar.gz | wc -l | tr -d ' ') entries"
cat dist/SHA256SUMS
