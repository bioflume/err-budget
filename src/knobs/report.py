"""Reports: the result as Markdown, and the run log as a table."""
from __future__ import annotations

import json
import os


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    return repr(v) if isinstance(v, str) else str(v)


def markdown(result) -> str:
    """The recommendation as Markdown (a settings table and the diagnostics)."""
    L = [f"## knobs: recommendation for `{result.observable}`", ""]
    if result.certified:
        L.append(f"**Target {result.target:.3g} met.** Estimated error "
                 f"{abs(result.error_estimate):.3g}.")
    else:
        L.append(f"**Target {result.target:.3g} not reachable with the values given.** "
                 f"Best estimate at the returned setting: {abs(result.error_estimate):.3g}.")
    L += ["", "| parameter | recommended | was |", "|---|---|---|"]
    for k, v in result.settings.items():
        was = result.start.get(k)
        L.append(f"| `{k}` | {_fmt(v)} | {'' if was == v else _fmt(was)} |")
    if result.limiting:
        L += ["", f"What limits the accuracy: `{result.limiting}`."]
    if result.budget_shares:
        shares = ", ".join(f"`{k}` {100 * s:.0f}%" for k, s in result.budget_shares.items())
        L += ["", f"Share of the error budget: {shares}."]
    if result.wasted:
        L += ["", "Refined past the point of having any effect:", ""]
        L += [f"- {s}" for s in result.wasted]
    if result.frozen:
        L += ["", "Frozen by the preliminary sweep (kept at your value): "
              + ", ".join(f"`{k}`" for k in result.frozen) + "."]
    if result.violations:
        L += ["", "Declarations contradicted by measurement:", ""]
        L += [f"- `{s}` was declared cost-only but moves the observable" for s in result.violations]
    ratio = result.calibration_cost / max(result.recommended_cost, 1e-300)
    L += ["", f"{result.runs_used} solver runs used; calibration cost {ratio:.0f}x one run "
          f"at the recommended setting."]
    if result.trace:
        L += ["", "<details><summary>Search path</summary>", "", "| step | setting |", "|---|---|"]
        for i, t in enumerate(result.trace):
            L.append(f"| {i} | " + ", ".join(f"{k}={_fmt(v)}" for k, v in t.items()) + " |")
        L += ["", "</details>"]
    return "\n".join(L) + "\n"


def run_table(log_path: str, limit: int | None = None) -> str:
    """The run log written by the command-line runner (runs.jsonl) as a text table."""
    from .runner import _read_log
    rows = _read_log(log_path)
    if limit:
        rows = rows[-limit:]
    if not rows:
        return "(no runs yet)"
    names = list(dict.fromkeys(n for r in rows for n in r.get("params", {})))
    head = f"{'run':<5} " + " ".join(f"{n:>12}" for n in names) + f" {'value':>14} {'cost':>10} {'status':<10}"
    out = [head, "-" * len(head)]
    for r in rows:
        vals = " ".join(f"{_fmt(r['params'].get(n, '-')):>12}" for n in names)
        v = r.get("value")
        vs = f"{v:>14.8g}" if isinstance(v, (int, float)) and v is not None else f"{'-':>14}"
        out.append(f"{r.get('run_id', ''):<5} {vals} {vs} {float(r.get('cost', 0) or 0):>10.4g} {r.get('status', ''):<10}")
    return "\n".join(out)
