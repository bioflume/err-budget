"""Revision 4: the mixed-precision GROMACS build with version 1.2.1.

1. Version 1.2.1 without the round-off probe and with it, at the targets of both
   revisions (1e-4, 3e-5, 1e-5, 3e-6, 1e-6), graded against the double-precision
   reference.  (The sentence of Section 7.3 that quoted "6 of 9 and 6 of 6 accepted
   starts at 1e-5 and 3e-6" came from exp44_gromacs_floor.py, which runs version 1.1.)
2. Why Table 3 counts the declines at 3e-5 and 1e-5 as "reachable": the number of
   mixed-precision settings whose energy happens to lie within the target of the
   double-precision reference, against the round-off spread of the build.
    python exp78_gromacs_floor_v121.py      # writes data/gromacs_floor_v121.json
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cases as C
import pooled as P

SPREAD6 = 3.0e-5          # relative spread of the six converged settings (exp44, sec82_floor.txt)


def main():
    gd, gs = C._gromacs()
    QD = gd.exact["potential"]
    out = dict(floor_probe_abs=gs.floor, floor_probe_rel=gs.floor / abs(QD), spread6_rel=SPREAD6, rows=[])
    print(f"probe floor {gs.floor:.4f} kJ/mol = {gs.floor/abs(QD):.2e} relative; declines below "
          f"{10*gs.floor/abs(QD):.2e}; six-setting spread {SPREAD6:.1e} relative")
    good = [i for i in gs.grid() if gs.usable(i)]
    for tg in [1e-4, 3e-5, 1e-5, 3e-6, 1e-6]:
        eps = tg * abs(QD)
        within = [i for i in good if abs(gs.tensor[i][0]["potential"] - QD) <= eps]
        dev_within = [abs(gs.tensor[i][0]["potential"] - QD) / abs(QD) for i in within]
        row = dict(target=tg, settings=len(good), within_target_by_chance=len(within),
                   largest_error_of_those=max(dev_within) if dev_within else None)
        case = C.Case("GMX-s", "GROMACS mixed", gs.tensor, gs.LAD, gs.HS, gs.KINDS, dict(gs.exact),
                      {"potential": [tg]}, floor=gs.floor)
        for cfg in ("v1.2.1", "v1.2.1-noprobe", "v1.1", "v1.1-noprobe"):
            s = P.summarise(P.replay(case, "potential", tg, P.CONFIGS[cfg]))
            row[cfg] = {k: s[k] for k in ("accepted", "missed", "declined_reach", "declined_unreach",
                                          "distinct", "worst", "pen_med", "reasons")}
        out["rows"].append(row)
        print(f"target {tg:.0e}: {len(within)}/{len(good)} settings lie within it by chance; "
              + "; ".join(f"{cfg}: accepted {row[cfg]['accepted']} missed {row[cfg]['missed']} "
                          f"(worst {row[cfg]['worst']:.2f} eps) distinct {row[cfg]['distinct']}"
                          for cfg in ("v1.2.1", "v1.2.1-noprobe", "v1.1", "v1.1-noprobe")), flush=True)
    json.dump(out, open(os.path.join(C.D, "gromacs_floor_v121.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
