"""Item: is F_s = 1.25 specific to L2?  Calibrate on each problem in turn.

Uses the margin scan of exp29 (data/margin_scan.json).  For each problem used
as the calibration set, F_s is the smallest scanned margin that is fully safe on
every one of its cells; the other problems are held out and graded at that F_s.
"""
import json
import os

import numpy as np

D = os.path.join(os.path.dirname(__file__), "data")
res = json.load(open(os.path.join(D, "margin_scan.json")))
FS = sorted(float(f) for f in res[0]["row"])
probs = sorted({r["cell"].split()[0] for r in res})
print(f"scanned margins: {FS}")
print(f"{'calibrate on':<14} {'F_s chosen':>10}   held-out cells: fully safe / all, "
      f"pooled safety, median of worst-case penalty")
for cal in probs:
    cells = [r for r in res if r["cell"].startswith(cal + " ")]
    need = [r["fmin"] for r in cells]
    if any(f is None for f in need):
        print(f"{cal:<14} no scanned margin is safe on every cell"); continue
    F = max(need)
    held = [r for r in res if not r["cell"].startswith(cal + " ")]
    ok = sum(r["fmin"] is not None and r["fmin"] <= F for r in held)
    row = [r["row"][str(F) if str(F) in r["row"] else f"{F:g}"] for r in held]
    cert = sum(x["certified"] for x in row); safe = sum(x["safe"] for x in row)
    pmax = np.median([x["pen_max"] for x in row])
    print(f"{cal:<14} {F:>10.2f}   {ok}/{len(held)}, {100*safe/max(cert,1):.2f}% "
          f"({cert-safe} misses in {cert} certified), {pmax:.1f}")
print("\nsmallest fully safe margin per cell:")
for r in res:
    print(f"  {r['cell']:<34} {r['fmin']}")
