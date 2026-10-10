"""Summary of pooled_<config>.json files: totals and pooled penalty statistics.

Aggregation (stated in every caption of the article): "median penalty" is the
median over cells of the per-cell median penalty; "geo-mean" the geometric mean
of the per-cell medians.
"""
import json, sys, os
import numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
for cfg in sys.argv[1:]:
    p = os.path.join(D, f"pooled_{cfg}.json")
    if not os.path.exists(p):
        print(cfg, "missing"); continue
    R = json.load(open(p))
    t = lambda k: sum(r.get(k, 0) for r in R)
    pm = [r["pen_med"] for r in R if r["accepted"] and np.isfinite(r["pen_med"])]
    cells_miss = sum(1 for r in R if r["missed"])
    usable = t('n') - t('start_failed')
    print(f"{cfg:<26} cells {len(R)} replays {t('n')} start-failed {t('start_failed')} "
          f"accepted {t('accepted')} ({100*t('accepted')/max(usable,1):.1f}% of usable starts) "
          f"missed {t('missed')} in {cells_miss} cells, declined-reachable {t('declined_reach')}, "
          f"declined-unreachable {t('declined_unreach')}, unconverged accepted {t('unconv')}, "
          f"low-order accepted {t('low_order')} (missed {t('low_order_missed')}), "
          f"median of cell medians {np.median(pm):.2f}, geo-mean {np.exp(np.mean(np.log(pm))):.2f}, "
          f"runs median {np.median([r['runs_med'] for r in R if np.isfinite(r['runs_med'])]):.0f}")
