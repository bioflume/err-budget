"""Check that the program frozen for the held-out tests is the one PREREGISTRATION.md fixed.

PREREGISTRATION.md ("Frozen code hash") fixes the program as the git tag v1.2.0-frozen at
commit 41d3db3, with the SHA-256 of the lines `<sha256>  <file name>` over src/knobs/*.py
in sorted order. This script recomputes that hash for

  1. paper/frozen_src/knobs/, the copy that the replays import, and
  2. src/knobs/ at the tag, and checks that the tag points at the registered commit.

    python tools/check_frozen.py                # both; fails if the tag is missing
    python tools/check_frozen.py --no-tag       # only the copy in paper/frozen_src
"""
import argparse
import hashlib
import io
import os
import re
import subprocess
import sys
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG = "v1.2.0-frozen"


def registered():
    text = open(os.path.join(ROOT, "PREREGISTRATION.md"), encoding="utf-8").read()
    section = text.split("## Frozen code hash", 1)[1]
    commit = re.search(r"\(commit ([0-9a-f]{7,40})\)", section).group(1)
    digest = re.search(r"\b([0-9a-f]{64})\b", section).group(1)
    return commit, digest


def tree_hash(files):
    """files: {name: bytes}; the hash of the sorted `<sha256>  <name>` lines."""
    lines = "".join(f"{hashlib.sha256(files[n]).hexdigest()}  {n}\n" for n in sorted(files))
    return hashlib.sha256(lines.encode()).hexdigest()


def frozen_copy():
    d = os.path.join(ROOT, "paper", "frozen_src", "knobs")
    return {n: open(os.path.join(d, n), "rb").read()
            for n in os.listdir(d) if n.endswith(".py")}


def tagged():
    git = ["git", "-C", ROOT]
    ok = subprocess.run(git + ["rev-parse", "-q", "--verify", f"refs/tags/{TAG}"],
                        capture_output=True, text=True)
    if ok.returncode:
        return None, None
    commit = subprocess.run(git + ["rev-list", "-n", "1", TAG], capture_output=True,
                            text=True, check=True).stdout.strip()
    tar = subprocess.run(git + ["archive", TAG, "src/knobs"], capture_output=True,
                         check=True).stdout
    files = {}
    with tarfile.open(fileobj=io.BytesIO(tar)) as tf:
        for m in tf.getmembers():
            if m.isfile() and os.path.dirname(m.name) == "src/knobs" and m.name.endswith(".py"):
                files[os.path.basename(m.name)] = tf.extractfile(m).read()
    return commit, files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-tag", action="store_true")
    args = ap.parse_args()
    commit, digest = registered()
    bad = []
    h = tree_hash(frozen_copy())
    print(f"paper/frozen_src/knobs: {h}")
    if h != digest:
        bad.append("paper/frozen_src/knobs does not match the registered hash")
    if not args.no_tag:
        tc, files = tagged()
        if tc is None:
            bad.append(f"the tag {TAG} is missing; the author pushes it with "
                       f"`git push origin {TAG}`")
        else:
            ht = tree_hash(files)
            print(f"{TAG} ({tc[:7]}) src/knobs: {ht}")
            if not tc.startswith(commit):
                bad.append(f"{TAG} points at {tc[:7]}, not at the registered commit {commit}")
            if ht != digest:
                bad.append(f"src/knobs at {TAG} does not match the registered hash")
    if bad:
        sys.exit("\n".join(bad))
    print(f"frozen program matches PREREGISTRATION.md ({digest[:12]}...)")


if __name__ == "__main__":
    main()
