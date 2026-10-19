"""Turn the margin scan into the LaTeX table, distinguishing the two ways a
cell can fail to deliver: certifying something that misses, and refusing
everything.  A scan that reports only safe/certified conflates them, and at
large F_s the method refuses rather than errs.
"""
import json
import os
import sys

D = os.path.join(os.path.dirname(__file__), "data")
rows = json.load(open(os.path.join(D, "margin_scan.json")))
FS = ["1.0", "1.05", "1.1", "1.15", "1.2", "1.25", "1.35", "1.5", "2.0", "3.0"]


def cellstr(r, f):
    d = r["row"].get(f) or r["row"].get(f + "0") or r["row"].get(f.rstrip("0"))
    if d is None:
        return "--"
    n, c, s = d["n"], d["certified"], d["safe"]
    if c == 0:
        return r"\emph{ref}"                 # refused every start
    if s == c == n:
        return r"$\checkmark$"               # certified all, all attained
    if s == c:
        return f"{100*c//n}\\%r"             # some refused, none missed
    return f"\\textbf{{{100*s//c}}}"         # certified some that missed


print(r"\begin{table}[t]\centering")
print(r"\caption{Safety against the margin $F_s$, with $F_s$ selected on the L2"
      r" cells and L3/L4 held out. $\checkmark$: every start certified and every"
      r" recommendation attained its target. \emph{ref}: the method refused every"
      r" start (a false refusal, not an error). A bold number is the percentage"
      r" of \emph{certified} recommendations that attained the target, i.e. a"
      r" genuine failure. $F_s^{\min}$ is the smallest margin on this scan with"
      r" no failure and no blanket refusal.}")
print(r"\label{tab:margin}")
print(r"\small")
print(r"\begin{tabular}{ll" + "r" * len(FS) + "r}")
print(r"\toprule")
print("cell & role & " + " & ".join("$" + f + "$" for f in FS)
      + r" & $F_s^{\min}$ \\")
print(r"\midrule")
lastrole = None
for r in rows:
    if lastrole is not None and r["role"] != lastrole:
        print(r"\midrule")
    lastrole = r["role"]
    cells = [cellstr(r, f) for f in FS]
    ok = []
    for f in FS:
        d = r["row"].get(f)
        if d and d["certified"] == d["n"] and d["safe"] == d["n"]:
            ok.append(float(f))
    fmin = f"${min(ok):.2f}$" if ok else "---"
    name = r["cell"].replace("_", r"\_").replace("rel", " rel")
    print(f"{name} & {r['role']} & " + " & ".join(cells) + f" & {fmin} " + r"\\")
print(r"\bottomrule")
print(r"\end{tabular}")
print(r"\end{table}")

cal = [r for r in rows if r["role"] == "calibration" and r["fmin"]]
held = [r for r in rows if r["role"] == "held out" and r["fmin"]]
if cal and held:
    print()
    print("%% calibration cells need at most F_s =", max(r["fmin"] for r in cal))
    print("%% held-out cells need at most F_s =", max(r["fmin"] for r in held))
