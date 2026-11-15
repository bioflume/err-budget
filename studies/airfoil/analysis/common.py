"""Shared helpers of the analyses: the run store, the configuration, the engine.

Every analysis reads the run JSONs of a store directory, either $AIRFOIL_WORK/runs on
Hamilton or the `runs/` of an unpacked tarball on the laptop, and the configuration
the runs were made with (`config/` of the code, or `packs/meta-*/config` of the
tarball), so the sweeps resolve to the same tags as on the cluster.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.dont_write_bytecode = True

from af import study  # noqa: E402
from af.cgrid import MESH_VERSION  # noqa: E402

ARTICLE_SRC = os.path.join(ROOT, "..", "..", "src")


def use_config(path):
    """Point the study definition at a configuration directory (None: keep)."""
    if path:
        study.CONFIG_DIR = os.path.abspath(path)
    return study.CONFIG_DIR


def find_unpacked(root):
    """(runs dir, config dir) of an unpacked tarball, or of $AIRFOIL_WORK."""
    runs = os.path.join(root, "runs")
    metas = sorted(glob.glob(os.path.join(root, "packs", "meta-*", "config")))
    return runs, (metas[-1] if metas else None)


class Store:
    """The run JSONs and error files of one store directory."""

    def __init__(self, runs_dir):
        self.dir = runs_dir
        self.results, self.errors, self.stale = {}, {}, []
        for p in glob.glob(os.path.join(runs_dir, "*.json")):
            name = os.path.basename(p)[:-5]
            if ".error-" in name:
                self.errors.setdefault(name.split(".error-")[0], []).append(p)
            elif not name.endswith(".tmp"):
                with open(p) as fh:
                    r = json.load(fh)
                # results made with an older mesh or set-up (PROTOCOL.md 14, A1-A3) are not results
                if r.get("mesh_version", 1) == MESH_VERSION and r.get("setup_version", 1) == study.SETUP_VERSION:
                    self.results[name] = r
                else:
                    self.stale.append(name)

    def get(self, run):
        return self.results.get(run.tag)

    def value(self, run, key="ld"):
        r = self.get(run)
        return math.nan if r is None else float(r.get(key, math.nan))

    def usable(self, run):
        """A converged result (the only kind any procedure may use)."""
        r = self.get(run)
        return bool(r and r.get("status") == "ok" and r.get("converged")
                    and math.isfinite(float(r.get("ld", math.nan))))


def reference_value(store, run, ff_run):
    """A5: the reference value of a design, L/D(d) + Delta_R, with Delta_R from the
    reference's first run (n_surf 91 at R = 5000, settled; in its result) and the far-field
    run (the same at R = 500).  Returns (value, info); value is nan unless all are usable."""
    info = dict(raw=math.nan, delta=math.nan, ld_5000=math.nan, ld_500=math.nan)
    if not store.usable(run):
        return math.nan, info
    j = store.get(run)
    info["raw"] = float(j["ld"])
    seq = j.get("sequence") or {}
    info["ld_5000"] = float(seq.get("ld", math.nan)) if seq.get("converged") else math.nan
    info["ld_500"] = store.value(ff_run) if store.usable(ff_run) else math.nan
    if not (math.isfinite(info["ld_5000"]) and math.isfinite(info["ld_500"])):
        return math.nan, info
    info["delta"] = study.far_field_correction(info["ld_5000"], info["ld_500"])
    return info["raw"] + info["delta"], info


def engine():
    """errbudget's engine: the vendored copy of the tarball, $ERRBUDGET_SRC, or the article."""
    for cand in (os.path.join(ROOT, "engine"), os.environ.get("ERRBUDGET_SRC", ""), ARTICLE_SRC):
        if cand and os.path.isdir(os.path.join(cand, "errbudget")):
            if cand not in sys.path:
                sys.path.insert(0, cand)
            break
    import errbudget
    return errbudget


def write_json(path, obj):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=_jsonable)


def _jsonable(x):
    try:
        import numpy as np
        if isinstance(x, np.bool_):
            return bool(x)
        if isinstance(x, np.integer):
            return int(x)
        if isinstance(x, np.floating):
            return float(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
    except ImportError:
        pass
    if isinstance(x, (set, tuple)):
        return list(x)
    raise TypeError(type(x))


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def fmt(x, nd=4):
    if x is None:
        return "-"
    if isinstance(x, float):
        if not math.isfinite(x):
            return "nan"
        return f"{x:.{nd}g}"
    return str(x)
