"""Revision 4: write RESULTS_SUMMARY.md and pack every result the rewrite needs.

    python exp82_collect_v6.py [archive.tar.gz]      # default ../knobs_v6_results.tar.gz

RESULTS_SUMMARY.md (in the repository root) records what was run and with which program,
and the headline numbers read back from the result files: the pooled table for every
configuration (with the refined and the revision-3 2D-2 study), the held-out replays, the
accounting of Table 2, the baselines, and the revision-4 analyses (ABINIT reference,
GROMACS floor, PPPM mechanism, relTol, unconverged counts).  A section that cannot be read
says so instead of stopping the script.  The archive holds the summary, every result file
of paper/data that revision 4 reads or writes, the regenerated tables and figures, the logs
of the analysis steps (logs/v6/) and hamilton/results.
"""
import glob
import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import tarfile
import time
from contextlib import redirect_stdout

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
D = os.path.join(HERE, "data")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

CONFIG_ORDER = ["v1.2.1", "v1.2", "v1.1", "v1.1+gated", "v1.1-noprobe", "v1.2.1-noprobe"]


def _json(name):
    p = os.path.join(D, name)
    return json.load(open(p)) if os.path.exists(p) else None


def versions():
    out = []
    try:
        import knobs
        src = os.path.dirname(knobs.__file__)
        sums = "".join(hashlib.sha256(open(f, "rb").read()).hexdigest() + "  " + os.path.basename(f) + "\n"
                       for f in sorted(glob.glob(os.path.join(src, "*.py"))))
        out.append(f"- knobs {getattr(knobs, '__version__', '?')} from {src} "
                   f"(sha256 of the sources {hashlib.sha256(sums.encode()).hexdigest()[:16]}...)")
    except Exception as exc:                       # noqa: BLE001
        out.append(f"- knobs: not importable ({exc})")
    try:
        g = subprocess.run(["git", "describe", "--tags", "--always", "--dirty"], cwd=ROOT,
                           capture_output=True, text=True)
        out.append(f"- git: {g.stdout.strip() or 'no git history'}")
    except OSError as exc:
        out.append(f"- git: not available ({exc})")
    import scipy
    out.append(f"- Python {platform.python_version()}, NumPy {np.__version__}, SciPy {scipy.__version__}, "
               f"{platform.node()}, {time.strftime('%Y-%m-%d %H:%M')}")
    return out


def pooled_totals():
    lines = ["| configuration | 2D-2 study | cells | replays | start failed | accepted | distinct | "
             "missed | declined (reachable) | declined (unreachable) | unconverged accepted | "
             "median penalty | geo. mean | median runs |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    files = sorted(glob.glob(os.path.join(D, "pooled_*.json")))
    names = [os.path.basename(f)[7:-5] for f in files]
    names = sorted(names, key=lambda n: (CONFIG_ORDER.index(n) if n in CONFIG_ORDER else 99, n))
    for name in names:
        R_all = json.load(open(os.path.join(D, f"pooled_{name}.json")))
        has_h = any(r["case"] == "2D-2h" for r in R_all)
        for study in (["refined", "rev3"] if has_h else ["rev3"]):
            drop = "2D-2" if study == "refined" else "2D-2h"
            R = [r for r in R_all if r["case"] != drop]
            s = lambda k: sum(r.get(k, 0) or 0 for r in R)
            pen = [r["pen_med"] for r in R if r.get("pen_med") is not None and np.isfinite(r["pen_med"])]
            runs = [r["runs_med"] for r in R if r.get("runs_med") is not None and np.isfinite(r["runs_med"])]
            lines.append(f"| {name} | {study} | {len(R)} | {s('n')} | {s('start_failed')} | {s('accepted')} | "
                         f"{s('distinct')} | {s('missed')} | {s('declined_reach')} | {s('declined_unreach')} | "
                         f"{s('unconv')} | {np.median(pen) if pen else float('nan'):.2f} | "
                         f"{np.exp(np.mean(np.log(pen))) if pen else float('nan'):.2f} | "
                         f"{np.median(runs) if runs else float('nan'):.0f} |")
    return lines


def misses(name="v1.2.1"):
    R = _json(f"pooled_{name}.json") or []
    rows = [r for r in R if r.get("missed")]
    if not rows:
        return [f"- {name}: no accepted setting missed in any cell."]
    return [f"- {name}: {r['case']} {r['ob']} {r['target']:g}: missed {r['missed']} of {r['accepted']} "
            f"(distinct missed {r.get('distinct_missed')}, worst {r['worst']:.2f} x target)" for r in rows]


def heldout():
    out = []
    for fn, lab in (("heldout_replay.json", "frozen v1.2 and v1.1"),
                    ("heldout_replay_current.json", "v1.2.1")):
        R = _json(fn)
        if R is None:
            out.append(f"- {fn}: missing")
            continue
        out.append(f"\n{lab} ({fn}):\n")
        out.append("| case | config | obs | target | starts | start failed | accepted | distinct | missed | "
                   "declined (reach.) | declined (unreach.) | penalty | runs | worst |")
        out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in R:
            out.append(f"| {r['case']} | {r['config']} | {r['ob']} | {r['target']:g} | {r['n']} | "
                       f"{r['start_failed']} | {r['accepted']} | {r.get('distinct', '')} | {r['missed']} | "
                       f"{r['declined_reach']} | {r['declined_unreach']} | {r['pen_med']:.2f} | "
                       f"{r['runs_med']:.0f} | {r['worst']:.2f} |")
    return out


def small(name, keys=None):
    m = _json(name)
    if m is None:
        return [f"- {name}: missing"]
    if keys and isinstance(m, dict):
        m = {k: m.get(k) for k in keys}
    txt = json.dumps(m, indent=1, default=str)
    return ["```json", txt[:6000] + ("\n... (truncated; full file in the archive)" if len(txt) > 6000 else ""),
            "```"]


def baselines():
    try:
        import exp72_baselines_pooled as B
    except Exception as exc:                       # noqa: BLE001
        return [f"- baselines: cannot import exp72 ({exc})"]
    doc = _json("baselines_v6.json")
    if not doc:
        return ["- baselines_v6.json: missing"]
    buf = io.StringIO()
    with redirect_stdout(buf):
        B.summary(doc["cells"])
    return [f"configuration of the program: {doc.get('config')}", "```", buf.getvalue().strip(), "```"]


def hamilton():
    res = os.path.join(ROOT, "hamilton", "results")
    out = []
    for j in ("j1", "j2", "j3", "j4", "j5"):
        p = os.path.join(res, f"{j}.txt")
        if os.path.exists(p):
            txt = open(p, errors="replace").read()
            out.append(f"- {j}.txt: {len(txt.splitlines())} lines (in the archive); last lines:")
            out.append("```")
            out.append("\n".join(txt.strip().splitlines()[-12:]))
            out.append("```")
        else:
            out.append(f"- {j}.txt: missing")
    return out


def main():
    arch = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "knobs_v6_results.tar.gz")
    S = ["# Results for revision 4 (knobs v1.2.1)", "", "## Program and machine", ""] + versions()
    sections = [
        ("Pooled table (Table 3) and ablations (Table 4): totals per configuration", pooled_totals),
        ("Misses of the default configuration", misses),
        ("Held-out problems (Table 5)", heldout),
        ("Accounting of the settings (Table 2)", lambda: small("accounting_v6.json")),
        ("ABINIT reference uncertainty (exp75)", lambda: small("abinit_reference_v6.json",
                                                               ["ref", "stored_ref", "reproduced", "unc",
                                                                "unc_exp45_rule", "parts"])),
        ("GROMACS round-off floor with v1.2.1 (exp78)", lambda: small("gromacs_floor_v121.json")),
        ("LAMMPS PPPM order and grid (exp74)", lambda: small("pppm_mechanism_v6.json")),
        ("2D-1 relTol and absolute tolerances (exp77; one row per residualControl x GAMG tolerance)",
         lambda: small("reltol_2d1.json")),
        ("Accepting unconverged runs, start by start (exp79)", lambda: small("unconverged_v6.json", ["total", "cases"])),
        ("Baselines (Table 6)", baselines),
        ("Hamilton jobs", hamilton),
    ]
    for title, fn in sections:
        S += ["", f"## {title}", ""]
        try:
            S += fn()
        except Exception as exc:                   # noqa: BLE001 - report and go on
            S.append(f"(could not be summarised: {type(exc).__name__}: {exc})")
    logs = sorted(glob.glob(os.path.join(ROOT, "logs", "v6", "*.txt")))
    S += ["", "## Logs of the analysis steps", ""] + [f"- logs/v6/{os.path.basename(f)}" for f in logs]
    summ = os.path.join(ROOT, "RESULTS_SUMMARY.md")
    open(summ, "w").write("\n".join(S) + "\n")
    print("written", summ)
    pats = ["RESULTS_SUMMARY.md", "paper/data/*_v6.json", "paper/data/pooled_*.json",
            "paper/data/heldout_replay*.json", "paper/data/gromacs_floor_v121.json",
            "paper/data/reltol_2d1.json", "paper/data/baselines_v6.json",
            "paper/figures/*.tex", "paper/figures/fig_search*", "paper/figures/*v6*",
            "logs/v6/*", "hamilton/results/*"]
    # the raw per-task files of the cluster (hamilton/results/raw.tar.gz, ~100 MB) stay behind
    files = sorted({f for p in pats for f in glob.glob(os.path.join(ROOT, p))
                    if os.path.basename(f) != "raw.tar.gz"})
    with tarfile.open(arch, "w:gz") as tf:
        for f in files:
            tf.add(f, arcname=os.path.relpath(f, ROOT))
    print(f"packed {len(files)} files into {arch} ({os.path.getsize(arch) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
