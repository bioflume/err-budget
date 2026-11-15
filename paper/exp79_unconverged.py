"""Revision 4 (referee point on the unconverged counts): what "accepting unconverged runs"
changes, start by start.

Table 4 reports that accepting unconverged runs adds about 1,500 acceptances, while the
text says that 39 of the accepted settings would then be unconverged runs.  Both are
true; they count different things.  For every admitted cell of every case that has
unconverged runs (the only cases where the two configurations can differ), this replays
version 1.2.1 and version 1.2.1 with unconverged runs accepted from every usable start
and tabulates, per case and in total:
  starts            usable starts (runs that produced a value)
  own run unconv.   starts whose own run did not converge (version 1.2.1 reports them
                    as "start failed" instead of searching from them)
  acc. 1.2.1        accepted by version 1.2.1
  acc. accept-unc.  accepted with unconverged runs accepted, split into acceptances from
                    starts whose own run did not converge and from converged starts
  unconv. accepted  accepted settings that are themselves unconverged runs (over starts,
                    distinct, and from converged starts)
  missed            misses with unconverged runs accepted
  changed           converged starts whose outcome differs between the two configurations

    python exp79_unconverged.py [CASE ...]
    ERRBUDGET_TABLE3_2D2=refined (default) | rev3, as for Table 3
Writes data/unconverged_v6.json.
"""
import json
import os
import sys

import cases as C
import pooled as P

A, B = "v1.2.1", "v1.2.1-accept-unconverged"


def per_start(case, ob, tg, cfg, a):
    o = P.replay(case, ob, tg, P.CONFIGS[cfg], starts=[a])
    if o["accepted"]:
        q = next(iter(o["proposals"]))
        return dict(acc=True, q=q, unconv=not case.converged(q), miss=bool(o["missed"]))
    return dict(acc=False, why=next(iter(o["reasons"]), None))


def main():
    keys = set(sys.argv[1:]) or None
    cases = C.all_cases()
    if any(c.key == "2D-2h" for c in cases):
        drop = "2D-2" if os.environ.get("ERRBUDGET_TABLE3_2D2", "refined") == "refined" else "2D-2h"
        cases = [c for c in cases if c.key != drop]
    out, tot = {}, {}
    hdr = (f"{'case':<7} {'cells':>5} {'starts':>7} {'own unc.':>8} {'acc 1.2.1':>9} {'acc a-u':>8} "
           f"{'from unc.':>9} {'from conv.':>10} {'unc.acc':>7} {'distinct':>8} {'unc.acc conv.st':>15} "
           f"{'missed':>6} {'changed':>7}")
    print(hdr)
    for case in cases:
        if keys and case.key not in keys:
            continue
        G = case.grid()
        if all(case.converged(i) for i in G if case.usable(i)):
            continue
        row = dict(cells=0, starts=0, own_unconv=0, acc_a=0, acc_b=0, acc_b_from_unconv=0,
                   acc_b_from_conv=0, unconv_acc=0, unconv_acc_from_conv=0, missed_b=0,
                   changed_conv=0, unconv_settings=set())
        starts = [i for i in G if case.usable(i)]
        for ob, tgs in case.targets.items():
            for tg in tgs:
                if not case.admitted(ob, tg):
                    continue
                row["cells"] += 1
                for a in starts:
                    ra, rb = per_start(case, ob, tg, A, a), per_start(case, ob, tg, B, a)
                    conv = case.converged(a)
                    row["starts"] += 1
                    row["own_unconv"] += not conv
                    row["acc_a"] += ra["acc"]
                    if rb["acc"]:
                        row["acc_b"] += 1
                        row["acc_b_from_unconv" if not conv else "acc_b_from_conv"] += 1
                        row["missed_b"] += rb["miss"]
                        if rb["unconv"]:
                            row["unconv_acc"] += 1
                            row["unconv_acc_from_conv"] += conv
                            row["unconv_settings"].add((ob, tg, rb["q"]))
                    if conv and (ra["acc"] != rb["acc"] or (ra["acc"] and ra["q"] != rb["q"])):
                        row["changed_conv"] += 1
        row["unconv_distinct"] = len(row.pop("unconv_settings"))
        out[case.key] = row
        for k, v in row.items():
            tot[k] = tot.get(k, 0) + v
        print(f"{case.key:<7} {row['cells']:>5} {row['starts']:>7} {row['own_unconv']:>8} {row['acc_a']:>9} "
              f"{row['acc_b']:>8} {row['acc_b_from_unconv']:>9} {row['acc_b_from_conv']:>10} "
              f"{row['unconv_acc']:>7} {row['unconv_distinct']:>8} {row['unconv_acc_from_conv']:>15} "
              f"{row['missed_b']:>6} {row['changed_conv']:>7}", flush=True)
    if tot:
        print(f"{'TOTAL':<7} {tot['cells']:>5} {tot['starts']:>7} {tot['own_unconv']:>8} {tot['acc_a']:>9} "
              f"{tot['acc_b']:>8} {tot['acc_b_from_unconv']:>9} {tot['acc_b_from_conv']:>10} "
              f"{tot['unconv_acc']:>7} {tot['unconv_distinct']:>8} {tot['unconv_acc_from_conv']:>15} "
              f"{tot['missed_b']:>6} {tot['changed_conv']:>7}")
        print(f"\nextra acceptances with unconverged runs accepted: {tot['acc_b'] - tot['acc_a']} = "
              f"{tot['acc_b_from_unconv']} from starts whose own run did not converge + "
              f"{tot['acc_b_from_conv'] - tot['acc_a']} from converged starts; of all "
              f"{tot['acc_b']} acceptances, {tot['unconv_acc']} ({tot['unconv_distinct']} distinct "
              f"settings) are unconverged runs, {tot['unconv_acc_from_conv']} of them from converged starts")
    json.dump(dict(configs=[A, B], cases=out, total=tot),
              open(os.path.join(C.D, "unconverged_v6.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
