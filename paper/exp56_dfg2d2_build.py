"""Build the OpenFOAM grid for the DFG benchmark 2D-2 (Re = 100, periodic).

Restartable and parallel: every setting is a separate task whose result is
written to its own JSON file, so any number of workers (on one machine, or as a
SLURM array) can run at once and an interrupted build resumes where it stopped.

    python exp56_dfg2d2_build.py list            # print the tasks
    python exp56_dfg2d2_build.py spin NA         # spin up mesh level NA
    python exp56_dfg2d2_build.py task I          # run task I (spins its mesh if needed)
    python exp56_dfg2d2_build.py worker          # run tasks until none is left
    python exp56_dfg2d2_build.py status          # progress and estimated time left
    python exp56_dfg2d2_build.py collect         # assemble data/dfg2d2_tensor.pkl

Grid: n_arc 5 x dt 4 x outer_tol 4 x p_tol 3 = 240 settings, plus 3 reference
runs on finer meshes and time steps (the reference and its uncertainty).
Work directory: $KNOBS_DFG2_WORK (default paper/work_dfg2d2).
"""
import fcntl
import itertools
import json
import os
import pickle
import socket
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import openfoam_dfg2d2 as F

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("KNOBS_DFG2_WORK", os.path.join(HERE, "work_dfg2d2"))
# results of revision 3 (corrected wall-cell centroids) go to results_v12; the
# spun-up meshes are shared, since only the post-processing changed
RES = os.path.join(WORK, os.environ.get("KNOBS_DFG2_RESULTS", "results_v12"))
NAMES = ["n_arc", "dt", "outer_tol", "p_tol"]
LAD = [[12, 16, 24, 32, 48], [0.004, 0.002, 0.001, 0.0005],
       [1e-2, 1e-3, 1e-4, 1e-5], [1e-4, 1e-6, 1e-8]]
GRID = list(itertools.product(*[range(len(l)) for l in LAD]))
REFS = [dict(n_arc=64, dt=0.00025, outer_tol=1e-5, p_tol=1e-8),     # the reference
        dict(n_arc=64, dt=0.0005, outer_tol=1e-5, p_tol=1e-8),      # its time-step error
        dict(n_arc=48, dt=0.00025, outer_tol=1e-5, p_tol=1e-8)]     # its mesh error
SPIN = dict(dt=0.001, outer_tol=1e-4, p_tol=1e-8)
CACHE = os.path.join(HERE, "data", os.environ.get("KNOBS_DFG2_CACHE", "dfg2d2_tensor_v12.pkl"))


def tasks():
    """Cheapest first, so that results accumulate early; references last."""
    def est(i):
        na, dt = LAD[0][i[0]], LAD[1][i[1]]
        return (na / 16.0) ** 2 / dt
    out = [("grid", i, {NAMES[k]: LAD[k][j] for k, j in enumerate(i)})
           for i in sorted(GRID, key=est)]
    out += [("ref", r, th) for r, th in enumerate(REFS)]
    return out


def _path(kind, key):
    tag = "_".join(map(str, key)) if kind == "grid" else f"ref{key}"
    return os.path.join(RES, f"{tag}.json")


class _Lock:
    """An exclusive, non-blocking lock on a file; released if the process dies."""
    def __init__(self, name):
        os.makedirs(os.path.join(WORK, "locks"), exist_ok=True)
        self.f = open(os.path.join(WORK, "locks", name), "a+")

    def acquire(self, block=False):
        try:
            fcntl.flock(self.f, fcntl.LOCK_EX | (0 if block else fcntl.LOCK_NB))
            self.f.seek(0); self.f.truncate()
            self.f.write(f"{socket.gethostname()} {os.getpid()} {time.ctime()}\n"); self.f.flush()
            return True
        except BlockingIOError:
            return False

    def release(self):
        fcntl.flock(self.f, fcntl.LOCK_UN)


def solver():
    return F.DFG2D2(workdir=os.path.join(WORK, "cases"), spin=SPIN)


def spun(na):
    return os.path.exists(os.path.join(WORK, "cases", f"mesh_{na}", "SPUN_UP"))


def spin(na, block=True):
    """Spin up mesh level na.  With block=False, return False at once if another
    worker is already spinning it (the caller then picks other work)."""
    if spun(na):
        return True
    lk = _Lock(f"spin_{na}")
    if not lk.acquire(block=block):
        return False
    try:
        t0 = time.time()
        solver().mesh(int(na))
        print(f"mesh {na} ready ({time.time()-t0:.0f}s)", flush=True)
        return True
    finally:
        lk.release()


def run_task(i):
    kind, key, th = tasks()[int(i)]
    out = _path(kind, key)
    if os.path.exists(out):
        return "done"
    if not spin(th["n_arc"], block=False):
        return "busy"                          # its mesh is being spun up elsewhere
    lk = _Lock(f"task_{i}")
    if not lk.acquire():
        return "busy"
    try:
        obs, cost, wall, info = solver().run(th)
        rec = dict(kind=kind, key=list(key) if kind == "grid" else key, theta=th, obs=obs,
                   cost=cost, wall=wall, info=info, host=socket.gethostname(),
                   finished=time.ctime())
        os.makedirs(RES, exist_ok=True)
        tmp = out + ".tmp"
        json.dump(rec, open(tmp, "w"))
        os.replace(tmp, out)
        print(f"task {i} {th}: {obs} cpu {cost:.0f}s wall {wall:.0f}s", flush=True)
        return "ran"
    finally:
        lk.release()


def worker():
    while True:
        ran = False
        for i, (kind, key, th) in enumerate(tasks()):
            if os.path.exists(_path(kind, key)):
                continue
            if run_task(i) == "ran":
                ran = True
                break
        if not ran:
            left = [i for i, (k, key, _) in enumerate(tasks()) if not os.path.exists(_path(k, key))]
            if not left:
                print("all tasks done", flush=True)
                return
            time.sleep(60)                     # the rest are held by other workers


def status():
    T = tasks()
    done = [json.load(open(_path(k, key))) for k, key, _ in T if os.path.exists(_path(k, key))]
    n = len(T)
    print(f"{len(done)}/{n} tasks done ({sum(1 for d in done if d['kind']=='grid')}/{len(GRID)} "
          f"grid, {sum(1 for d in done if d['kind']=='ref')}/{len(REFS)} reference)")
    spun = [na for na in LAD[0] + [64] if os.path.exists(
        os.path.join(WORK, "cases", f"mesh_{na}", "SPUN_UP"))]
    print(f"meshes spun up: {spun}")
    if done:
        # crude estimate: wall time per task scales like cells / dt
        rate = [d["wall"] / ((d["theta"]["n_arc"] / 16.0) ** 2 / d["theta"]["dt"]) for d in done]
        r = sorted(rate)[len(rate) // 2]
        left = sum((th["n_arc"] / 16.0) ** 2 / th["dt"] for k, key, th in T
                   if not os.path.exists(_path(k, key)))
        print(f"estimated work left: {r*left/3600:.1f} worker-hours")
    bad = [d for d in done if not (d["cost"] < float("inf"))]
    print(f"failed runs among those done: {len(bad)}")


def collect():
    T = tasks()
    tensor, ref = {}, {}
    for kind, key, th in T:
        p = _path(kind, key)
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        info = {k: v for k, v in d["info"].items() if not k.endswith("series")}
        if kind == "grid":
            tensor[tuple(key)] = (d["obs"], d["cost"], info)
        else:
            ref[key] = (d["theta"], d["obs"])
    pickle.dump(dict(tensor=tensor, ref=ref, LAD=LAD, NAMES=NAMES), open(CACHE, "wb"))
    print(f"{len(tensor)}/{len(GRID)} grid settings and {len(ref)}/{len(REFS)} references -> {CACHE}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "worker"
    if cmd == "list":
        for i, (k, key, th) in enumerate(tasks()):
            print(i, k, th, "done" if os.path.exists(_path(k, key)) else "")
    elif cmd == "spin":
        spin(int(sys.argv[2]), block=True)
    elif cmd == "task":
        print(run_task(int(sys.argv[2])))
    elif cmd == "worker":
        worker()
    elif cmd == "status":
        status()
    elif cmd == "collect":
        collect()
    else:
        raise SystemExit(__doc__)
