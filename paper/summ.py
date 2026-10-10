"""Summary of pooled_<config>.json files: totals and pooled penalty statistics."""
import json, sys, os
import numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
for cfg in sys.argv[1:]:
    p = os.path.join(D, f"pooled_{cfg}.json")
    if not os.path.exists(p):
        print(cfg, "missing"); continue
    R = json.load(open(p))
    t = lambda k: sum(r[k] for r in R)
    pm = [r["pen_med"] for r in R if r["accepted"] and np.isfinite(r["pen_med"])]
    cells_miss = sum(1 for r in R if r["missed"])
    print(f"{cfg:<26} cells {len(R)} replays {t('n')} accepted {t('accepted')} "
          f"({100*t('accepted')/t('n'):.1f}%) missed {t('missed')} in {cells_miss} cells, "
          f"declined-reachable {t('declined_reach')}, unconverged accepted {t('unconv')}, "
          f"median of cell median penalties {np.median(pm):.2f}, "
          f"geo-mean {np.exp(np.mean(np.log(pm))):.2f}, runs median {np.median([r['runs_med'] for r in R]):.0f}")
