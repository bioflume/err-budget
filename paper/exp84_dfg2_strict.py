"""OpenFOAM 2D-2 (refined study, case 2D-2h): the stricter rule after a failed finer run
(S-5.2, "A stricter rule"; Section 8).

The released program charges a parameter whose finer run failed its last coarser change
(LocalBudget.per_knob_error).  Under the stricter rule a discretization parameter keeps that
charge only if the current value and the two coarser ones ran; otherwise its error cannot
be estimated and the program declines (unmeasured_axis).  Iterative-solver parameters keep
the released rule.  The rule is applied here by wrapping per_knob_error; the engine is not
changed.  Every test case of 2D-2h is searched from every completed setting with the
stricter rule and compared with the released program: the stored replay
(data/pooled_v1.2.1.json) for the graded test cases, a fresh one for the two ungraded
targets of S-5.2.

    python exp84_dfg2_strict.py        # about two minutes on four cores; writes data/dfg2_strict.json
"""
import json
import os
from multiprocessing import Pool

import numpy as np

import cases as C
import pooled as P
from errbudget import LocalBudget
from paths import DATA, OUT

UNGRADED = (("cd_max", 3e-3), ("cl_max", 1e-2))
_per_knob_error = LocalBudget.per_knob_error


def strict_per_knob_error(self, p_grid=None):
    _per_knob_error(self, p_grid)
    for k in range(len(self.m)):
        a = self.anchor[k]
        if not (self._disc(k) and self.fine_failed[k]) or k in self.unmeasured:
            continue
        if a < 2 or not (np.isfinite(self.e[k][a - 1]) and np.isfinite(self.e[k][a - 2])):
            self.unmeasured.append(k)
            self.rk[k] = np.inf
    self.E_knob = float(np.nansum(self.rk))
    return self.E_knob


def run(job):
    ob, tg, strict = job
    LocalBudget.per_knob_error = strict_per_knob_error if strict else _per_knob_error
    o = P.replay(C._dfg2_refined(), ob, tg, P.CONFIGS["v1.2.1"])
    return ob, tg, strict, dict(accepted=o["accepted"], missed=o["missed"],
                                declined=o["declined_reach"] + o["declined_unreach"])


if __name__ == "__main__":
    case = C._dfg2_refined()
    tests = [(ob, tg) for ob, tgs in case.targets.items() for tg in tgs
             if case.admitted(ob, tg) or (ob, tg) in UNGRADED]
    jobs = [(ob, tg, True) for ob, tg in tests] + [(ob, tg, False) for ob, tg in UNGRADED]
    with Pool(int(os.environ.get("ERRBUDGET_NPROC", 4))) as pool:
        res = {(ob, tg, s): r for ob, tg, s, r in pool.map(run, jobs)}
    for row in json.load(open(os.path.join(DATA, "pooled_v1.2.1.json"))):
        if row["case"] == "2D-2h":
            res[(row["ob"], row["target"], False)] = dict(
                accepted=row["accepted"], missed=row["missed"],
                declined=row["declined_reach"] + row["declined_unreach"])
    print("2D-2h, released rule against the stricter rule (accepted / missed / declined)")
    out = []
    for ob, tg in tests:
        a, b = res[(ob, tg, False)], res[(ob, tg, True)]
        graded = case.admitted(ob, tg)
        print(f"  {ob:<7} {tg:<6g} {'graded  ' if graded else 'ungraded'}  released "
              f"{a['accepted']:5d} / {a['missed']:3d} / {a['declined']:5d}   stricter "
              f"{b['accepted']:5d} / {b['missed']:3d} / {b['declined']:5d}   lost "
              f"{a['accepted'] - b['accepted']}")
        out.append(dict(ob=ob, target=tg, graded=graded, released=a, stricter=b))
    print(f"misses removed: {sum(r['released']['missed'] - r['stricter']['missed'] for r in out)}")
    json.dump(out, open(os.path.join(OUT, "dfg2_strict.json"), "w"), indent=1)
