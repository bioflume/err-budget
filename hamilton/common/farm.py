"""A restartable task farm for SLURM arrays (no locks, no database).

A job module defines

    TASKS()          -> list of dicts, each with a unique "id" (str)
    run(task, out)   -> runs one task, writes its result into the directory `out`
                        (a JSON file named result.json at least)

and this script runs them:

    python -m farm <module> list                 # print tasks and their status
    python -m farm <module> array K NK [key=value ...]
                                                 # array element K of NK: tasks i = K, K+NK, ...
                                                 # of those whose fields match the filters
    python -m farm <module> one ID               # one task by id
    python -m farm <module> status               # counts

Results go to $WORK/<module>/<id>/.  A task whose result.json exists is skipped,
so arrays can be resubmitted after time-outs; a task that raised writes
error.txt and is retried on the next submission.
"""
import importlib
import json
import os
import sys
import time
import traceback


def _dir(mod, tid):
    return os.path.join(os.environ.get("WORK", "."), mod.__name__.split(".")[-1], tid)


def done(mod, t):
    return os.path.exists(os.path.join(_dir(mod, t["id"]), "result.json"))


def run_one(mod, t):
    out = _dir(mod, t["id"])
    os.makedirs(out, exist_ok=True)
    if done(mod, t):
        return "done"
    err = os.path.join(out, "error.txt")
    if os.path.exists(err):
        # from an earlier attempt: kept as error_prev.txt, so that a later time-out is not
        # reported with an old error
        os.replace(err, os.path.join(out, "error_prev.txt"))
    with open(os.path.join(out, "attempts.txt"), "a") as fa:
        fa.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} job {os.environ.get('SLURM_JOB_ID', '-')} "
                 f"array {os.environ.get('SLURM_ARRAY_JOB_ID', '-')}_{os.environ.get('SLURM_ARRAY_TASK_ID', '-')}\n")
    t0 = time.time()
    try:
        mod.run(t, out)
        if os.path.exists(os.path.join(out, "error.txt")):
            os.remove(os.path.join(out, "error.txt"))
        print(f"[{time.strftime('%H:%M:%S')}] {t['id']}: finished in {time.time()-t0:.0f}s",
              flush=True)
        return "ran"
    except Exception:
        open(os.path.join(out, "error.txt"), "w").write(traceback.format_exc())
        print(f"[{time.strftime('%H:%M:%S')}] {t['id']}: FAILED (see {out}/error.txt)",
              flush=True)
        return "failed"


def main():
    modname, cmd = sys.argv[1], sys.argv[2]
    mod = importlib.import_module(modname)
    T = mod.TASKS()
    if cmd == "list":
        for i, t in enumerate(T):
            print(i, t["id"], "done" if done(mod, t) else "")
    elif cmd == "status":
        n = sum(done(mod, t) for t in T)
        err = sum(os.path.exists(os.path.join(_dir(mod, t["id"]), "error.txt")) for t in T)
        started = sum(not done(mod, t) and not os.path.exists(os.path.join(_dir(mod, t["id"]), "error.txt"))
                      and os.path.exists(os.path.join(_dir(mod, t["id"]), "attempts.txt")) for t in T)
        print(f"{modname}: {n}/{len(T)} done, {err} with errors, {started} started and unfinished "
              f"(running, or stopped by a time-out)")
    elif cmd == "array":
        k, nk = int(sys.argv[3]), int(sys.argv[4])
        for f in sys.argv[5:]:
            key, val = f.split("=", 1)
            T = [t for t in T if str(t.get(key)) == val]
        for i in range(k, len(T), nk):
            run_one(mod, T[i])
    elif cmd == "count":
        for f in sys.argv[3:]:
            key, val = f.split("=", 1)
            T = [t for t in T if str(t.get(key)) == val]
        print(len(T))
    elif cmd == "one":
        t = [t for t in T if t["id"] == sys.argv[3]][0]
        print(run_one(mod, t))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
