"""Revision 4: how the grid sizes of Table 2 become the start counts of Table 3.

For every case: the number of settings in the exhaustive grid (the product of the
numbers of candidate values, which is printed); how many runs produced no value
(crashed, diverged, out of memory: never a start, never proposed); how many produced
a value but did not converge (counted as starts whose own run failed, "start
failed"); and how many converged.  Also the candidate values of every parameter, so
that Table 2 can list m_1 x ... x m_p.
    python exp76_accounting.py            # writes data/accounting_v6.json
"""
import json
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C


def reason(info):
    """Why a run produced no value, from the message its driver stored."""
    u = str(info.get("unavailable", "")) if isinstance(info, dict) else ""
    if not u:
        return "no value (no message stored)"
    if "GAMGSolver.C" in u:
        return "GAMG could not be set up (FOAM FATAL ERROR in GAMGSolver)"
    if "sigFpe" in u or "Floating point" in u or "printStack" in u:
        return "floating-point exception in the solver (diverged)"
    if "at ??:?" in u:      # a stack trace whose first lines were cut off when it was stored
        return "solver crashed (stack trace only; signal not stored)"
    if "out of memory" in u.lower() or "MemoryError" in u:
        return "out of memory"
    if "timeout" in u.lower():
        return "time limit"
    if "non-finite" in u:
        return "non-finite result"
    return u.strip().splitlines()[0][:60]


def account(case):
    G = case.grid()
    nval = [i for i in G if case.usable(i)]
    novalue = [i for i in G if not case.usable(i)]
    unconv = [i for i in nval if not case.converged(i)]
    why = Counter()
    for i in novalue:
        info = case.tensor[i][2] if len(case.tensor[i]) > 2 else {}
        why[reason(info)] += 1
    return dict(case=case.key, label=case.label, dims=[len(l) for l in case.LAD],
                settings=len(G), produced_value=len(nval), no_value=len(novalue),
                not_converged=len(unconv), converged=len(nval) - len(unconv),
                starts=len(nval), start_failed=len(unconv),
                ladders=[[float(x) if isinstance(x, (int, float, np.floating, np.integer)) else str(x)
                          for x in l] for l in case.LAD],
                kinds=list(case.KINDS), models=list(case.models) if case.models else None,
                reference=dict(case.exact), uncertainty=dict(case.unc),
                admitted={ob: [tg for tg in tgs if case.admitted(ob, tg)] for ob, tgs in case.targets.items()},
                not_admitted={ob: [tg for tg in tgs if not case.admitted(ob, tg)] for ob, tgs in case.targets.items()},
                no_value_reasons=dict(why.most_common(5)), cost=case.cost_note)


def main():
    rows = []
    loaders = [C._manufactured, lambda: [C._dfg1()], lambda: [C._dfg2()], lambda: [C._dfg2_refined()],
               C._lammps, C._gromacs, lambda: [C._abinit()], lambda: [C._qe()]]
    for L in loaders:
        try:
            got = L()
        except FileNotFoundError as exc:
            print(f"skipping: {exc}")
            continue
        for case in got:
            a = account(case)
            rows.append(a)
            print(f"{a['case']:<7} {' x '.join(map(str, a['dims'])):<22} = {a['settings']:>5} settings; "
                  f"no value {a['no_value']:>4}; value but not converged {a['not_converged']:>4} "
                  f"(start failed); converged {a['converged']:>5}; starts {a['starts']:>5}; "
                  f"admitted {a['admitted']}; not admitted {a['not_admitted']}", flush=True)
            if a["no_value"]:
                print(f"        no-value reasons: {a['no_value_reasons']}")
    json.dump(rows, open(os.path.join(C.D, "accounting_v6.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
