## Arm 1: PRE sees exactly errbudget's runs (608 searches per case)

| Test case | errbudget missed / declined | PRE missed / declined | errbudget cost ratio | PRE cost ratio | errbudget K | PRE K | PRE coverage of Q* (±2s) | median |m−Q*|/ε |
|---|---|---|---|---|---|---|---|---|
| c_D 0.01 | 0 / 0 | 52 / 0 | 7.7 | 3.2 | 113 | 1,365 | 30% | 0.19 |
| c_D 0.003 | 0 / 0 | 1 / 559 | 8.8 | 1 | 91.8 | 1,613 | 8% | 0.74 |
| c_D 0.001 | 0 / 0 | 12 / 549 | 6.6 | 1 | 76.5 | 540 | 7% | 2.2 |
| Δp 0.03 | 0 / 0 | 0 / 0 | 3 | 1 | 444 | 1,414 | 15% | 0.058 |
| Δp 0.01 | 0 / 0 | 4 / 3 | 7.9 | 3.2 | 224 | 549 | 16% | 0.23 |
| c_L 0.1 | 0 / 0 | 0 / 276 | 1.6 | 1 | 2,227 | 559 | 90% | 0.11 |
| c_L 0.03 | 0 / 0 | 0 / 416 | 14.2 | 1.6 | 385 | 776 | 70% | 0.26 |

PRE at z = 1 and z = 3 (missed / declined, median cost ratio):

| Test case | z = 1 | z = 3 | coverage ±1s / ±3s | 2s ≤ ε | |m−Q*| ≤ ε |
|---|---|---|---|---|---|
| c_D 0.01 | 82 / 0, 1.6 | 17 / 0, 3.2 | 11% / 46% | 100% | 100% |
| c_D 0.003 | 16 / 398, 18 | 1 / 561, 1 | 1% / 14% | 100% | 99% |
| c_D 0.001 | 18 / 541, 1 | 12 / 566, 2.8 | 0% / 9% | 63% | 10% |
| Δp 0.03 | 0 / 0, 1 | 0 / 0, 1 | 8% / 23% | 100% | 100% |
| Δp 0.01 | 4 / 0, 3.2 | 1 / 45, 3.2 | 3% / 30% | 100% | 100% |
| c_L 0.1 | 0 / 224, 1 | 0 / 295, 1 | 85% / 96% | 55% | 82% |
| c_L 0.03 | 0 / 281, 1.6 | 0 / 448, 1.6 | 63% / 88% | 32% | 81% |

## Arm 2: PRE's own design, budget B x1 errbudget's median search cost (50 seeds per case)

| Test case | budget (cells×iter.) | errbudget missed / declined | PRE missed / declined | errbudget cost ratio | PRE cost ratio | errbudget-3p cost ratio | errbudget K | PRE K | PRE coverage (±2s) | median |m−Q*|/ε |
|---|---|---|---|---|---|---|---|---|---|---|
| c_D 0.01 | 2.95e+08 | 0/608 / 0 | 1/50 / 0 | 7.7 | 1 | 28.3 | 113 | 2,723 | 48% | 0.096 |
| c_D 0.003 | 8.1e+08 | 0/608 / 0 | 6/50 / 2 | 8.8 | 1 | — | 91.8 | 995 | 76% | 0.23 |
| c_D 0.001 | 2.54e+09 | 0/608 / 0 | 0/50 / 0 | 6.6 | 2.8 | — | 76.5 | 179 | 100% | 0.47 |
| Δp 0.03 | 5.19e+07 | 0/608 / 0 | 0/50 / 5 | 3 | 1 | 3 | 444 | 1,447 | 20% | 0.074 |
| Δp 0.01 | 1.78e+08 | 0/608 / 0 | 0/50 / 2 | 7.9 | 3.2 | 7.9 | 224 | 543 | 30% | 0.19 |
| c_L 0.1 | 1.96e+08 | 0/608 / 0 | 0/50 / 42 | 1.6 | 1 | 1.6 | 2,227 | 5,504 | 90% | 0.14 |
| c_L 0.03 | 3.84e+08 | 0/608 / 0 | 0/50 / 49 | 14.2 | 1 | 48.7 | 385 | 11,325 | 96% | 0.44 |

## Arm 2: PRE's own design, budget B x3 errbudget's median search cost (50 seeds per case)

| Test case | budget (cells×iter.) | errbudget missed / declined | PRE missed / declined | errbudget cost ratio | PRE cost ratio | errbudget-3p cost ratio | errbudget K | PRE K | PRE coverage (±2s) | median |m−Q*|/ε |
|---|---|---|---|---|---|---|---|---|---|---|
| c_D 0.01 | 8.85e+08 | 0/608 / 0 | 0/50 / 0 | 7.7 | 1 | 28.3 | 113 | 8,641 | 86% | 0.061 |
| c_D 0.003 | 2.43e+09 | 0/608 / 0 | 0/50 / 0 | 8.8 | 1 | — | 91.8 | 3,038 | 100% | 0.16 |
| c_D 0.001 | 7.61e+09 | 0/608 / 0 | 0/50 / 0 | 6.6 | 2.8 | — | 76.5 | 200 | 100% | 0.46 |
| Δp 0.03 | 1.56e+08 | 0/608 / 0 | 0/50 / 1 | 3 | 1 | 3 | 444 | 4,330 | 38% | 0.062 |
| Δp 0.01 | 5.35e+08 | 0/608 / 0 | 0/50 / 0 | 7.9 | 3.2 | 7.9 | 224 | 1,638 | 40% | 0.17 |
| c_L 0.1 | 5.87e+08 | 0/608 / 0 | 0/50 / 49 | 1.6 | 1 | 1.6 | 2,227 | 16,619 | 98% | 0.13 |
| c_L 0.03 | 1.15e+09 | 0/608 / 0 | 0/50 / 50 | 14.2 | — | 48.7 | 385 | — | 88% | 11 |

## Headline criteria (protocol section 6, z = 2)

```
{
 "arm1": {
  "criterion1_less_reliable": true,
  "pre_misses_where_errbudget_does_not": [
   "cd 0.01",
   "cd 0.001"
  ],
  "coverage_below_90": [
   "cd 0.01",
   "cd 0.003",
   "cd 0.001",
   "dp 0.03",
   "dp 0.01",
   "cl 0.03"
  ],
  "criterion2_cheaper_same_reliability": false,
  "cases_cost_ratio_1p5x_lower": [
   "cd 0.01",
   "cd 0.003",
   "cd 0.001",
   "dp 0.03",
   "dp 0.01",
   "cl 0.1",
   "cl 0.03"
  ],
  "criterion3_less_usable": false,
  "cases_declining_over_10pct": [
   "cd 0.003",
   "cd 0.001",
   "cl 0.1",
   "cl 0.03"
  ]
 },
 "arm2": {
  "criterion1_less_reliable": true,
  "pre_misses_where_errbudget_does_not": [
   "cd 0.01",
   "cd 0.003"
  ],
  "coverage_below_90": [
   "cd 0.01",
   "cd 0.003",
   "dp 0.03",
   "dp 0.01"
  ],
  "criterion2_cheaper_same_reliability": false,
  "cases_cost_ratio_1p5x_lower": [
   "cd 0.01",
   "cd 0.003",
   "cd 0.001",
   "dp 0.03",
   "dp 0.01",
   "cl 0.1",
   "cl 0.03"
  ],
  "criterion3_less_usable": false,
  "cases_declining_over_10pct": [
   "cl 0.1",
   "cl 0.03"
  ]
 }
}
```
