#!/usr/bin/env python3
"""Verdict of the 15-minute check job (RUNBOOK.md, step 4).

    python3 analysis/check_smoke.py --work $AIRFOIL_WORK

Passes when both tiny runs of the `check` sweep converged with finite forces, y+ was
read from the wall, the mesh has the cells the generator promised, the solver is
v2606 and echoed the decay control of kOmegaSST, and the airFoil2D tutorial file has
the tolerances that setting (b) of PROTOCOL.md 7 cites.  Prints CHECK PASSED or CHECK
FAILED with the reasons.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import MESH_VERSION, Store, study  # noqa: E402

TUTORIAL = dict(p_tol=1e-6, p_reltol=0.1, u_tol=1e-8, res=1e-5)


def tutorial_values(path):
    """Tolerances of the airFoil2D fvSolution.  v2606 has separate U and nuTilda blocks;
    a combined "(U|nuTilda)" block is accepted too."""
    txt = open(path).read()
    solvers = txt[txt.index("solvers"):]

    def block(key):
        m = re.search(r'(?:^|\s)%s\s*\{([^}]*)\}' % key, solvers)
        return m.group(1) if m else None

    p = block("p")
    u = block("U") or block(r'"\(U\|nuTilda\)"')
    nut = block("nuTilda") or u
    res = re.search(r"residualControl\s*\{([^}]*)\}", txt).group(1)
    num = lambda blk, k: float(re.search(r"\b%s\s+([-+0-9.eE]+?);" % k, blk).group(1))
    return dict(p_tol=num(p, "tolerance"), p_reltol=num(p, "relTol"),
                u_tol=max(num(u, "tolerance"), num(nut, "tolerance")),
                res=max(num(res, "p"), num(res, "U"), num(res, "nuTilda")))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.environ.get("AIRFOIL_WORK", "."))
    a = ap.parse_args(argv)
    store = Store(os.path.join(a.work, "runs"))
    problems = []
    # the mesh checks (PROTOCOL.md 14, A2)
    import json
    for run in study.sweep_meshcheck():
        p = os.path.join(a.work, "meshcheck", f"{run.tag}.json")
        if not os.path.exists(p):
            problems.append(f"mesh check {run.tag}: not run")
            continue
        m = json.load(open(p))
        if m.get("mesh_version") != MESH_VERSION:
            problems.append(f"mesh check {run.tag}: made with an older mesh version")
        cm = m.get("checkMesh", {})
        print(f"mesh {run.tag}: {m.get('cells')} cells, checkMesh non-orthogonality {cm.get('max_non_orth')} "
              f"(generator {m.get('check', {}).get('max_non_orth_deg', float('nan')):.1f}), skewness "
              f"{cm.get('max_skewness')}, aspect {cm.get('max_aspect')}: {'passed' if m.get('passed') else 'FAILED'}")
        if not m.get("passed"):
            problems.append(f"mesh check {run.tag}: {cm.get('failures') or m.get('error')}")
    for i, run in enumerate(study.sweep_check()):
        r = store.get(run)
        if r is None:
            errs = store.errors.get(run.tag, [])
            problems.append(f"check {i}: no result" + (f" ({len(errs)} error files: {errs[-1]})" if errs else ""))
            continue
        name = f"check {i} ({run.tag})"
        if r["status"] != "ok" or not r["converged"]:
            problems.append(f"{name}: status {r['status']}, stop {r['stop_reason']}")
        for k in ("cl", "cd", "ld"):
            if not math.isfinite(float(r.get(k, math.nan))):
                problems.append(f"{name}: {k} is not a number (forceCoeffs output?)")
        if "yplus" not in r:
            problems.append(f"{name}: no y+ ({r.get('errors')})")
        cm = r["mesh"]["checkMesh"]
        if cm.get("cells") != r["mesh"]["generator"]["cells"]:
            problems.append(f"{name}: checkMesh {cm.get('cells')} cells, generator {r['mesh']['generator']['cells']}")
        ver = str(r["openfoam"].get("version"))
        if "2606" not in ver:
            problems.append(f"{name}: OpenFOAM version {ver}, not v2606")
        if not r["openfoam"].get("decay_control_echoed"):
            problems.append(f"{name}: kOmegaSST did not echo decayControl (free-stream decay not controlled)")
        print(f"{name}: L/D {r['ld']:.5g}, cl {r['cl']:.5g}, cd {r['cd']:.5g}, {r['iterations']} iterations, "
              f"{r['cells']} cells, cpu {r['cost']['cpu_s']:.1f} s, y+ mean {r.get('yplus', {}).get('mean', math.nan):.3g} "
              f"max {r.get('yplus', {}).get('max', math.nan):.3g}, checkMesh: max non-orth "
              f"{cm.get('max_non_orth')}, {cm.get('failed_checks')} failed checks")
    tut = os.path.join(a.work, "tutorial_airFoil2D_fvSolution")
    if os.path.exists(tut):
        try:
            v = tutorial_values(tut)
            print(f"airFoil2D tutorial: {v}")
            for k, want in TUTORIAL.items():
                if not math.isclose(v[k], want, rel_tol=1e-9):
                    problems.append(f"tutorial {k} = {v[k]}, PROTOCOL.md 7(b) cites {want}")
        except Exception as exc:
            problems.append(f"tutorial file not understood: {exc}")
    else:
        problems.append(f"tutorial file missing ({tut}): run slurm/setup.sh")
    if problems:
        print("CHECK FAILED:\n  " + "\n  ".join(problems))
        return 1
    print("CHECK PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
