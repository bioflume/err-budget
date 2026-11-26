"""Check that the human-owned files of this repository are unchanged.

The files listed in PROTECTED (one glob per line) are the requirements, the acceptance
tests and the record of the results. Their owner, not a coding agent, decides when they
change (AGENTS.md). PROTECTED.sha256 holds the SHA-256 of every file the globs match;
CI runs this script and fails if any file was changed, removed or added.

    python tools/check_protected.py            # check; exit 1 on any difference
    python tools/check_protected.py --update   # rewrite PROTECTED.sha256 (author only;
                                               # stage new protected files first)
"""
import argparse
import glob
import hashlib
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIST = os.path.join(ROOT, "PROTECTED")
MANIFEST = os.path.join(ROOT, "PROTECTED.sha256")


def tracked():
    """Files in the git index (committed or staged), so that ignored outputs of a local
    replay (figures in PDF, logs) are not mistaken for added protected files. Stage a new
    protected file before `--update`. Outside git: None."""
    r = subprocess.run(["git", "-C", ROOT, "ls-files", "-z"], capture_output=True)
    return set(r.stdout.decode().split("\0")) - {""} if r.returncode == 0 else None


def protected_files():
    files = set()
    with open(LIST) as fh:
        for line in fh:
            pattern = line.split("#", 1)[0].strip()
            if not pattern:
                continue
            matches = glob.glob(os.path.join(ROOT, pattern), recursive=True)
            if not matches:
                sys.exit(f"PROTECTED: the pattern {pattern!r} matches no file")
            files.update(os.path.relpath(m, ROOT) for m in matches if os.path.isfile(m))
    under_git = tracked()
    return sorted(f for f in files if os.path.basename(f) != ".DS_Store"
                  and "__pycache__" not in f.split(os.sep)
                  and (under_git is None or f in under_git))


def sha256(path):
    with open(os.path.join(ROOT, path), "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--update", action="store_true")
    args = ap.parse_args()
    current = {f: sha256(f) for f in protected_files()}
    if args.update:
        with open(MANIFEST, "w") as fh:
            for f, h in current.items():
                fh.write(f"{h}  {f}\n")
        print(f"wrote {len(current)} entries to PROTECTED.sha256")
        return
    stored = {}
    with open(MANIFEST) as fh:
        for line in fh:
            h, f = line.rstrip("\n").split("  ", 1)
            stored[f] = h
    changed = sorted(f for f in current if f in stored and current[f] != stored[f])
    added = sorted(set(current) - set(stored))
    removed = sorted(set(stored) - set(current))
    for label, fs in (("changed", changed), ("added", added), ("removed", removed)):
        for f in fs:
            print(f"{label}: {f}")
    if changed or added or removed:
        sys.exit("protected files differ from PROTECTED.sha256: their owner must approve "
                 "the change and run `python tools/check_protected.py --update`")
    print(f"{len(current)} protected files unchanged")


if __name__ == "__main__":
    main()
