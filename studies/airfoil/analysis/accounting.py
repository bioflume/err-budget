#!/usr/bin/env python3
"""Account for every run of a sweep (RUNBOOK.md, step 8 and the checks of each stage).

    python3 analysis/accounting.py --runs DIR [--config DIR] [--sacct FILE] SWEEP [SWEEP ...]

For each sweep: runs expected; results (converged; unconverged: stalled or at the
cap; failed: diverged); infrastructure errors (with their message) and runs given up
after three; runs without any record (missing or still running); the Slurm array
list to resubmit; and, with sacct's records, the job states (FAILED, TIMEOUT,
OUT_OF_MEMORY, CANCELLED) of the stage.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, study, use_config  # noqa: E402


def compress(idx):
    out, idx, i = [], sorted(idx), 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
            j += 1
        out.append(f"{idx[i]}" if i == j else f"{idx[i]}-{idx[j]}")
        i = j + 1
    return ",".join(out)


def account(store: Store, sweep: str) -> dict:
    runs = study.sweep(sweep)
    out = dict(sweep=sweep, expected=len(runs), converged=[], stalled=[], max_iter=[], other_unconverged=[],
               failed=[], errors={}, gave_up=[], missing=[])
    for i, r in enumerate(runs):
        res = store.get(r)
        if res is not None:
            if res.get("status") == "failed":
                out["failed"].append(i)
            elif res.get("converged"):
                out["converged"].append(i)
            elif res.get("stop_reason") == "stalled":
                out["stalled"].append(i)
            elif res.get("stop_reason") == "max_iter":
                out["max_iter"].append(i)
            else:
                out["other_unconverged"].append(i)
            continue
        errs = store.errors.get(r.tag, [])
        if len(errs) >= 3:
            out["gave_up"].append(i)
        elif errs:
            msgs = []
            for p in sorted(errs):
                with open(p) as fh:
                    e = json.load(fh)
                msgs.append(f"{e.get('stage')}: {e.get('message')}")
            out["errors"][i] = msgs
        else:
            out["missing"].append(i)
    out["results"] = sum(len(out[k]) for k in ("converged", "stalled", "max_iter", "other_unconverged", "failed"))
    out["resubmit"] = compress(out["missing"] + list(out["errors"]))
    out["complete"] = out["results"] + len(out["gave_up"]) == out["expected"]
    return out


def sacct_states(path):
    """Counts of job states of the array tasks in a `sacct -P` dump (steps excluded)."""
    if not path or not os.path.exists(path):
        return {}
    lines = open(path).read().strip().splitlines()
    if not lines:
        return {}
    hdr = lines[0].split("|")
    if "State" not in hdr or "JobID" not in hdr:
        return {}
    js, jid = hdr.index("State"), hdr.index("JobID")
    counts = {}
    for line in lines[1:]:
        f = line.split("|")
        if len(f) <= max(js, jid) or "." in f[jid]:
            continue
        s = f[js].split()[0] if f[js] else "?"
        counts[s] = counts.get(s, 0) + 1
    return counts


def report(a: dict) -> str:
    L = [f"{a['sweep']}: {a['expected']} runs expected, {a['results']} results "
         f"({len(a['converged'])} converged, {len(a['stalled'])} stalled, {len(a['max_iter'])} at the cap, "
         f"{len(a['other_unconverged'])} other unconverged, {len(a['failed'])} failed), "
         f"{len(a['errors'])} with errors, {len(a['gave_up'])} given up, {len(a['missing'])} without a record"]
    for i, msgs in sorted(a["errors"].items()):
        L.append(f"   error at index {i}: {msgs[-1]}")
    if a["gave_up"]:
        L.append(f"   given up (report these): {compress(a['gave_up'])}")
    if a["resubmit"]:
        L.append(f"   to resubmit: bash slurm/submit.sh {a['sweep']}   (indices {a['resubmit']})")
    L.append("   COMPLETE" if a["complete"] else "   NOT COMPLETE")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweeps", nargs="+")
    ap.add_argument("--runs", required=True)
    ap.add_argument("--config")
    ap.add_argument("--sacct")
    a = ap.parse_args(argv)
    use_config(a.config)
    store = Store(a.runs)
    for s in a.sweeps:
        print(report(account(store, s)))
    st = sacct_states(a.sacct)
    if st:
        print("sacct states of the jobs: " + ", ".join(f"{k} {v}" for k, v in sorted(st.items())))
        bad = {k: v for k, v in st.items() if k in ("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "OUT_OF_ME+",
                                                       "CANCELLED", "NODE_FAIL")}
        if bad:
            print("   look at the Slurm logs of these tasks (logs/af-<sweep>_<job>_<index>.out)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
