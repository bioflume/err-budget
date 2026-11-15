#!/usr/bin/env python3
"""The settling measurement (study.sweep_settle), on Hamilton or the laptop:

    python3 analysis/settle.py              # $AIRFOIL_WORK/runs, or --runs DIR

For each run (the reference mesh of m 4 %, t 12 % at R = 50,000 and 5,000, to 40,000
iterations without a residual stop): L/D at 10,000-iteration marks, and the iteration from
which L/D stays settled: from which on the relative spread of L/D over the preceding
`--window` iterations, (max - min) / L/D, never exceeds the tolerance again.  Tolerances
1e-3, 3e-4, 1e-4 and 3e-5; the question of the measurement is 1e-4 over 2,000."""
from __future__ import annotations

import argparse
import collections
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Store, fmt, md_table, study  # noqa: E402

TOLS = (1e-3, 3e-4, 1e-4, 3e-5)


def rolling_spread(x, w):
    """(max - min) over the w + 1 values ending at each index (monotone deques, O(n))."""
    out = np.full(len(x), np.nan)
    qmax, qmin = collections.deque(), collections.deque()
    for i, v in enumerate(x):
        while qmax and x[qmax[-1]] <= v:
            qmax.pop()
        qmax.append(i)
        while qmin and x[qmin[-1]] >= v:
            qmin.pop()
        qmin.append(i)
        while qmax[0] < i - w:
            qmax.popleft()
        while qmin[0] < i - w:
            qmin.popleft()
        if i >= w:
            out[i] = x[qmax[0]] - x[qmin[0]]
    return out


def settled_from(it, ld, window, tol):
    """First iteration from which the relative spread over `window` stays below tol, or None."""
    spread = rolling_spread(ld, window) / np.abs(ld)
    bad = np.nonzero(~(spread <= tol))[0]           # nan (too early) counts as not settled
    if len(bad) == len(ld):
        return None
    last_bad = bad[-1] if len(bad) else -1
    return None if last_bad == len(ld) - 1 else int(it[last_bad + 1])


def load(runs_dir, j):
    path = os.path.join(runs_dir, j.get("history_file") or os.path.join("history", f"{j['tag']}.csv"))
    a = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return a[:, 0], a[:, 1], a[:, 2]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(os.environ.get("AIRFOIL_WORK", "."), "runs"))
    ap.add_argument("--window", type=int, default=2000)
    a = ap.parse_args(argv)
    store = Store(a.runs)
    rows, marks = [], []
    for r in study.sweep_settle():
        j = store.get(r)
        label = f"R = {r.radius:g}"
        if j is None:
            errs = store.errors.get(r.tag, [])
            rows.append([label, "no result" + (f" ({len(errs)} errors)" if errs else "")] + [""] * (5 + len(TOLS)))
            continue
        it, cl, cd = load(a.runs, j)
        ld = cl / cd
        res = j.get("residuals_final") or {}
        seq = j.get("sequence") or {}
        settled = [settled_from(it, ld, a.window, tol) for tol in TOLS]
        rows.append([label, f"{j['status']}/{j['stop_reason']}", int(it[-1]), fmt(float(ld[-1]), 7),
                     fmt(max(res.values()), 2) if res else "", f"{seq.get('stop', '-')}, {seq.get('iterations', '-')}",
                     fmt(j["cost"]["core_s"] / 3600, 3)]
                    + [fmt(s) if s is not None else "not settled" for s in settled])
        marks.append([label] + [fmt(float(ld[min(len(ld) - 1, int(np.searchsorted(it, k)))]), 7)
                                for k in (10000, 20000, 30000, 40000)])
    print(f"settling of L/D: spread over the preceding {a.window} iterations, relative")
    print(md_table(["run", "outcome", "iterations", "L/D at the end", "largest final residual",
                    "first run (A4)", "core-hours"] + [f"settled to {t:g} from" for t in TOLS], rows))
    print()
    print(md_table(["run", "L/D at 10,000", "at 20,000", "at 30,000", "at 40,000"], marks))


if __name__ == "__main__":
    main()
