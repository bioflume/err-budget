#!/usr/bin/env python3
"""Run an errbudget job file through errbudget's Python API (RUNBOOK.md, step 10).

    python3 live_tune.py errbudget_job.yaml [--out DIR]

The article's repository has the engine and its Python API but not the `errbudget`
command line that README.md describes, so this runner reads the job file itself:
`command` with {{job_dir}} and {{parameter}} placeholders, `parameters` (name, kind,
values coarse to fine, start), `observable` (name, regex on stdout), `target`,
`relative`, `cost` (regex) and `failure.unconverged_if_regex`.  A non-zero exit or a
missing observable is a failed run (errbudget.RunFailed); a match of the failure
pattern is an unconverged run, which version 1.2.1 treats as failed.  It writes
errbudget_result.json and errbudget_result.txt to --out.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "analysis"))
from common import engine  # noqa: E402


def load(path):
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


def fill(command, job_dir, values):
    out = command.replace("{{job_dir}}", job_dir)
    for k, v in values.items():
        out = out.replace("{{%s}}" % k, repr(v) if isinstance(v, float) else str(v))
    if "{{" in out:
        raise ValueError(f"unfilled placeholder in: {out}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("job")
    ap.add_argument("--out", default=".")
    a = ap.parse_args(argv)
    eb = engine()
    job = load(a.job)
    job_dir = os.path.dirname(os.path.abspath(a.job))
    from common import study
    if os.path.basename(a.job) == "errbudget_job.yaml":
        lists = {p["name"]: [float(v) for v in p["values"]] for p in job["parameters"]}
        want = {n: [float(v) for v in vals] for n, vals in zip(study.NAMES, study.ladders_now())}
        if lists != want:
            sys.exit(f"the job file's lists {lists} differ from the registered ones {want} "
                     f"(after the far-field amendment of PROTOCOL.md 3.2, multiply the radius values by "
                     f"{(study.load_config('stage1.json') or {}).get('radius_factor', 1)})")
    params = [eb.Parameter(p["name"], p["values"], kind=p["kind"]) for p in job["parameters"]]
    start = {p["name"]: p["start"] for p in job["parameters"] if "start" in p}
    name = job["observable"]["name"]
    obs_re = re.compile(job["observable"]["regex"])
    cost_re = re.compile(job["cost"]["regex"]) if job.get("cost") else None
    fail_re = (re.compile(job["failure"]["unconverged_if_regex"])
               if job.get("failure", {}).get("unconverged_if_regex") else None)
    log = []

    def solver(**kw):
        cmd = fill(job["command"], job_dir, kw)
        t0 = time.time()
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        rec = dict(setting=kw, returncode=r.returncode, seconds=time.time() - t0, stdout=r.stdout[-400:])
        log.append(rec)
        if r.returncode != 0:
            raise eb.RunFailed(f"exit {r.returncode}: {r.stdout[-200:]} {r.stderr[-200:]}")
        m = obs_re.search(r.stdout)
        if not m:
            raise eb.RunFailed("observable not found in the output")
        val = float(m.group(1))
        cost = float(cost_re.search(r.stdout).group(1)) if cost_re else None
        conv = not (fail_re and fail_re.search(r.stdout))
        rec.update(value=val, cost=cost, converged=conv)
        out = {name: val, "converged": conv}
        return (out, cost) if cost is not None else out

    res = eb.tune(solver, params, target=float(job["target"]), observable=name, start=start,
                  relative=bool(job.get("relative", False)), verbose=True)
    print(res)
    os.makedirs(a.out, exist_ok=True)
    doc = dict(job=os.path.abspath(a.job), certified=res.certified, settings=res.settings, start=res.start,
               error_estimate=res.error_estimate, target=res.target, limiting=res.limiting,
               budget_shares=res.budget_shares, wasted=res.wasted, runs_used=res.runs_used,
               calibration_cost=res.calibration_cost, recommended_cost=res.recommended_cost,
               engine=getattr(eb, "__version__", None), runs=log)
    with open(os.path.join(a.out, "errbudget_result.json"), "w") as fh:
        json.dump(doc, fh, indent=1, default=float)
    with open(os.path.join(a.out, "errbudget_result.txt"), "w") as fh:
        fh.write(str(res) + "\n")
    print(f"written {os.path.join(a.out, 'errbudget_result.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
