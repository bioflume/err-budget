"""Held-out tests H1 and H2 (and H3 if its Hamilton results are present), replayed
with the frozen program exactly as fixed in PREREGISTRATION.md.

    python exp62_heldout.py            # the frozen program (git tag v1.2.0-frozen), as
                                       # registered, and version 1.1 for comparison
    python exp62_heldout.py --current  # the released program (version 1.2.1, the
                                       # amendment registered after H2 and before H3)

The frozen source is extracted from the git tag and its hash checked against
PREREGISTRATION.md before anything is replayed.

The program was called knobs when the held-out tests were registered, and the frozen
source keeps that name (src/knobs at the tag, paper/frozen_src/knobs).  The replay scripts
import the current name, errbudget, so in the frozen mode the frozen package is registered
under that name before any of them is imported.
"""
import glob
import hashlib
import json
import os
import pickle
import sys

import subprocess
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CURRENT = "--current" in sys.argv
if CURRENT:
    SRC = os.path.join(HERE, "..", "src", "errbudget")
else:
    # the frozen source from the git tag; in a release without the git history, the
    # copy in paper/frozen_src (its hash is checked below either way)
    tmp = tempfile.mkdtemp(prefix="errbudget_frozen_")
    arch = subprocess.run(["git", "archive", "v1.2.0-frozen", "src/knobs"], cwd=os.path.join(HERE, ".."),
                          capture_output=True)
    if arch.returncode == 0 and arch.stdout:
        subprocess.run(["tar", "-x", "-C", tmp], input=arch.stdout, check=True)
        sys.path.insert(0, os.path.join(tmp, "src"))
        SRC = os.path.join(tmp, "src", "knobs")
    else:
        sys.path.insert(0, os.path.join(HERE, "frozen_src"))
        SRC = os.path.join(HERE, "frozen_src", "knobs")
if CURRENT:
    import errbudget as program
else:
    import importlib
    program = importlib.import_module("knobs")               # the frozen package, by its name
    sys.modules["errbudget"] = program
    for sub in ("_core", "_asymptotic", "api", "problems"):
        mod = importlib.import_module(f"knobs.{sub}")
        sys.modules[f"errbudget.{sub}"] = mod
        setattr(program, sub, mod)
assert os.path.dirname(program.__file__) == os.path.abspath(SRC) or CURRENT, program.__file__
import cases as C
import pooled as P


def source_hash():
    sums = "".join(hashlib.sha256(open(f, "rb").read()).hexdigest() + "  " + os.path.basename(f) + "\n"
                   for f in sorted(glob.glob(os.path.join(SRC, "*.py"))))
    return hashlib.sha256(sums.encode()).hexdigest()


reg = open(os.path.join(HERE, "..", "PREREGISTRATION.md")).read()
h = source_hash()
if not CURRENT:
    if h not in reg:
        raise SystemExit(f"source hash {h} is not the frozen one recorded in PREREGISTRATION.md")
    print(f"frozen program verified (sha256 {h[:16]}...), from {SRC}")
else:
    print(f"current program, errbudget {program.__version__} (sha256 {h[:16]}...)")
CFGS = ("v1.2.1",) if CURRENT else ("v1.2", "v1.1")


def heldout_case(key):
    if key == "H1":
        import heldout_fem as M
        exact, tg = M.EXACT, {"mean": [1e-3, 3e-4], "flux": [3e-2, 1e-2, 3e-3]}
        label = "H1: 3D Poisson, scikit-fem"
    else:
        import heldout_cavity as M
        exact, tg = M.REF, {o: [3e-2, 1e-2, 3e-3] for o in M.REF}
        label = "H2: cavity Re 1000, OpenFOAM"
    t = pickle.load(open(os.path.join(C.D, f"heldout_{key}.pkl"), "rb"))
    t = {i: (v[0], v[1], dict(v[2])) for i, v in t.items()}
    return C.Case(key, label, t, M.LAD, M.HS, M.KINDS, dict(exact), tg)


rows = []
cases = [heldout_case(k) for k in ("H1", "H2")
         if os.path.exists(os.path.join(C.D, f"heldout_{k}.pkl"))]
j5 = os.path.join(HERE, "..", "hamilton", "results", "j5_tensor.pkl")
if os.path.exists(j5):
    sys.path.insert(0, os.path.join(HERE, "..", "hamilton", "J5_abc3d"))
    sys.path.insert(0, os.path.join(HERE, "..", "hamilton", "common"))
    import abc_hpc as A
    t = pickle.load(open(j5, "rb"))
    cases.append(C.Case("H3", "H3: 3D ABC flow, OpenFOAM (Hamilton)", t, A.LAD, A.HS, A.KINDS,
                        dict(A.EXACT), {"ek": [1e-3, 1e-4], "eps": [3e-2, 1e-2, 3e-3],
                                        "umax": [1e-3, 3e-4]}))
for case in cases:
    G = case.grid()
    print(f"\n{case.label}: {len(G)} settings, {sum(case.usable(i) for i in G)} usable, "
          f"{sum(case.usable(i) and not case.converged(i) for i in G)} not converged")
    fin = [len(l) - 1 for l in case.LAD]
    for ob in case.targets:
        e = []
        for a in range(len(case.LAD[0])):
            i = tuple([a] + fin[1:])
            e.append(abs(case.tensor[i][0][ob] - case.exact[ob]) / abs(case.exact[ob])
                     if case.usable(i) else np.nan)
        print(f"  {ob:<6} error against the first parameter at the finest others: "
              + "  ".join(f"{x:.2e}" for x in e))
    for cfg in CFGS:
        for ob, tgs in case.targets.items():
            for tg in tgs:
                s = P.summarise(P.replay(case, ob, tg, P.CONFIGS[cfg]))
                s.update(config=cfg, case=case.key)
                rows.append(s)
                print(f"  {cfg:<5} {ob:<6} {tg:<7.0e} starts {s['n']:<4} start-failed "
                      f"{s['start_failed']:<4} accepted {s['accepted']:<4} missed {s['missed']:<3} "
                      f"declined(reachable) {s['declined_reach']:<4} declined(unreachable) "
                      f"{s['declined_unreach']:<4} penalty {s['pen_med']:.2f}/{s['pen_max']:.2f} "
                      f"runs {s['runs_med']:.0f} worst {s['worst']:.2f}", flush=True)
json.dump(rows, open(os.path.join(C.D, "heldout_replay_current.json" if CURRENT
                             else "heldout_replay.json"), "w"), indent=1)
