"""The `knobs` command.

    knobs init [FILE]            write an annotated job file (default knobs.yaml)
    knobs check FILE             validate the job file and run the starting setting once
    knobs run FILE               tune; print the report; write FILE's output .json and .md
    knobs report RESULT.json     print a saved result again
    knobs runs FILE              show the run log of a job (what ran, what it gave)
    knobs version
"""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__, ENGINE_VERSION
from .api import TuneResult
from .jobfile import load_job, run_job, check_job, init_job
from .runner import JobError
from .report import run_table


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="knobs",
                                 description="Choose a solver's numerical parameters for a target "
                                             "accuracy in one observable, at close to minimal cost.")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("init", help="write an annotated job file to start from")
    p.add_argument("file", nargs="?", default="knobs.yaml")
    p.add_argument("--force", action="store_true", help="overwrite an existing file")

    p = sub.add_parser("check", help="validate a job file and run the starting setting once")
    p.add_argument("file")
    p.add_argument("--no-run", action="store_true", help="validate only, run nothing")

    p = sub.add_parser("run", help="tune the parameters described in a job file")
    p.add_argument("file")
    p.add_argument("--fresh", action="store_true",
                   help="ignore the run log and run everything again")
    p.add_argument("--quiet", action="store_true", help="do not print each run")
    p.add_argument("--no-write", action="store_true", help="print only; write no result files")

    p = sub.add_parser("report", help="print a saved result (knobs_result.json) again")
    p.add_argument("file")
    p.add_argument("--markdown", action="store_true")

    p = sub.add_parser("runs", help="show the run log of a job")
    p.add_argument("file", help="the job file (its workdir holds runs.jsonl)")
    p.add_argument("--last", type=int, default=None, help="only the last N runs")

    sub.add_parser("version", help="print the package and engine versions")

    a = ap.parse_args(argv)
    if a.cmd is None:
        ap.print_help()
        return 2
    try:
        return _dispatch(a)
    except JobError as exc:
        print(f"knobs: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nknobs: interrupted; `knobs run` again resumes from the run log", file=sys.stderr)
        return 130


def _dispatch(a) -> int:
    if a.cmd == "version":
        print(f"knobs {__version__} (engine {ENGINE_VERSION}, the article's version)")
        return 0
    if a.cmd == "init":
        path = init_job(a.file, overwrite=a.force)
        print(f"wrote {path}")
        print("next: edit it (the comments say what each part is), then")
        print(f"      knobs check {path}   and   knobs run {path}")
        return 0
    if a.cmd == "check":
        job = load_job(a.file)
        ok = check_job(job, run=not a.no_run)
        return 0 if ok else 1
    if a.cmd == "run":
        job = load_job(a.file)
        log = (lambda s: None) if a.quiet else print
        print(job.describe())
        print("")
        result = run_job(job, log=log, fresh=a.fresh, write=not a.no_write)
        print("")
        print(result)
        if not a.no_write:
            print(f"\nwritten: {job.output}.json and {job.output}.md"
                  f"   (run log: {os.path.join(job.workdir, 'runs.jsonl')})")
        return 0 if result.certified else 3
    if a.cmd == "report":
        r = TuneResult.load(a.file)
        print(r.markdown() if a.markdown else r)
        return 0
    if a.cmd == "runs":
        job = load_job(a.file)
        print(run_table(os.path.join(job.workdir, "runs.jsonl"), limit=a.last))
        return 0
    return 2


if __name__ == "__main__":                                   # python -m knobs.cli
    sys.exit(main())
