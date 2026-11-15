"""analysis/reuse.py end to end on a synthetic study: the article's stored 2D-1 runs at
Re = 20 (if the article repository is present), with the error scaled up with Re."""
import json
import os
import pickle
import sys

import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "analysis"))
from cyl import study  # noqa: E402

ARTICLE = os.path.join(HERE, "..", "..")
TENSOR = os.path.join(ARTICLE, "paper", "data", "dfg_tensor_v12.pkl")
pytestmark = pytest.mark.skipif(not os.path.exists(TENSOR), reason="article repository not present")
BENCH = 5.57953523384


def level(re_):
    return BENCH * (20.0 / re_) ** 0.6


def write_store(d):
    T = pickle.load(open(TENSOR, "rb"))
    os.makedirs(d, exist_ok=True)
    for re_ in study.MEMBERS:
        g = re_ / 20.0                                   # harder at larger Re
        for idx in study.GRID:
            r = study.run_at(re_, study.setting_of(idx))
            o, cost, info = T[idx]
            failed = "unavailable" in info
            rel = (o["cd"] - BENCH) / BENCH if not failed else 0.0
            x = dict(r.to_dict(), status="failed" if failed else "ok",
                     converged=bool(info.get("converged", False)) and not failed,
                     cd=float("nan") if failed else level(re_) * (1 + g * rel),
                     cpu=1.0 if failed else float(cost), iterations=info.get("iters") or 2000,
                     cells_iters=(info.get("cells") or 768) * (info.get("iters") or 2000))
            json.dump(x, open(os.path.join(d, f"{r.tag}.json"), "w"))
        for r in study.reference_runs(re_):
            e = 3.56e-4 * (64.0 / r.n_arc) ** 2 * g
            x = dict(r.to_dict(), status="ok", converged=True, cd=level(re_) * (1 + e), cpu=1000.0,
                     iterations=5000, cells_iters=1e9)
            json.dump(x, open(os.path.join(d, f"{r.tag}.json"), "w"))


def test_analysis_on_a_synthetic_study(tmp_path, capsys):
    import reuse
    runs = tmp_path / "runs"
    write_store(str(runs))
    assert reuse.main(["--runs", str(runs), "--out", str(tmp_path / "res"), "--starts-every", "40",
                       "--engine", os.path.join(ARTICLE, "src")]) == 0
    s = json.load(open(tmp_path / "res" / "reuse_summary.json"))
    assert s["accounting"]["grid"]["expected"] == 10920 and s["accounting"]["grid"]["missing"] == 0
    assert s["hardest"] == 40.0
    assert s["graded"] == {"0.01": True, "0.003": True}
    t1 = s["targets"]["0.01"]
    assert t1["a"]["certified"]
    met = [v["met"] for v in t1["reuse"]["(a)"].values()]
    assert len(met) == 13 and all(met)                    # easier members inherit the setting
    pb = t1["payback"]["cpu"]
    assert pb["K"] > 1 and pb["mu"] > 1 and pb["N_star"] > 0
    assert t1["every_start"]["starts"] > 0
    out = capsys.readouterr().out
    assert "hardest member: Re = 40" in out and "payback (cpu)" in out
