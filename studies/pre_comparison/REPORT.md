# errbudget against probabilistic Richardson extrapolation on OpenFOAM 2D-1

**Result.** Used to select a setting, PRE misses targets that errbudget meets from the same runs. Its 95% band
contains the exact value in only 7–30% of fits for drag and Δp. When it returns a setting, the setting is
1.6–8.6 times cheaper (median). The protocol was committed before grading (e1d1bce), and its criterion 1
("PRE less reliable") holds in both arms.

## PRE and its validation

- **Code and method.** The authors' code (NewcastleRSE/sparse-probabilistic-richardson-extrapolation,
  3918107) implements sparse PRE (SPRE) and the JRSS-B method (GRE). I used SPRE with their default
  white-noise kernel; the authors show GRE needs too many runs at d = 3.
- **Implementation.** Their code takes 1–10 minutes per fit, so `spre_white.py` restates the model in closed
  form.
- **Validation** (`results/validation.json`):
  - V0: their 15 unit tests against their MATLAB code pass.
  - V1: a rerun of their cubature experiment reproduces the published Fig. 1 estimates to 10 digits and the
    Fig. 2 errors to 3×10⁻⁴. The measured slopes are 2.00 and 3.99–4.00, the stated O(h²) and O(h⁴).
  - V2: the closed form matches their code on 9 of 10 data sets: same A, and m and s within 2×10⁻⁸ s. The
    tenth diverges above order 4 in an ill-conditioned basis, with m differing by 0.18 s. Capping the order
    at 4 leaves the findings unchanged (38 misses instead of 52 at c_D 10⁻²).

**Deviations.**
- The exploratory GRE check was cut to 3 data sets with a 15-minute cap, for the CPU limit. No fit finished,
  so GRE is untested here.
- Before the commit, V2 was cut from 21 data sets to 10 and Matérn-1/2 was dropped, both for CPU time.

**errbudget's replay**, graded by the same code as PRE, reproduces Table S-3.1 exactly and Table 7's K
(median 224).

## Arm 1: PRE sees exactly errbudget's runs (608 searches per case, z = 2)

| Case | Missed (eb / PRE) | Declined (eb / PRE) | Cost ratio (eb / PRE) | K (eb / PRE) | PRE coverage (m ± 2s) |
|---|---|---|---|---|---|
| c_D 10⁻² | 0 / **52** | 0 / 0 | 7.7 / 3.2 | 113 / 1,365 | 30% |
| c_D 3×10⁻³ | 0 / 1 | 0 / 559 | 8.8 / 1.0 | 92 / 1,613 | 8% |
| c_D 10⁻³ | 0 / **12** | 0 / 549 | 6.6 / 1.0 | 77 / 540 | 7% |
| Δp 3×10⁻² | 0 / 0 | 0 / 0 | 3.0 / 1.05 | 444 / 1,414 | 15% |
| Δp 10⁻² | 0 / 4 | 0 / 3 | 7.9 / 3.2 | 224 / 549 | 16% |
| c_L 10⁻¹ | 0 / 0 | 0 / 276 | 1.6 / 1.0 | 2,227 / 559 | 90% |
| c_L 3×10⁻² | 0 / 0 | 0 / 416 | 14.2 / 1.6 | 385 / 776 | 70% |

At z = 3, PRE still misses 17 searches at c_D 10⁻² and 12 at 10⁻³.

## Arm 2: PRE's own sequential design at errbudget's median search cost (50 designs per case, z = 2)

| Case | Missed (eb / PRE) | Declined (eb / PRE) | Cost ratio (eb / eb-3p / PRE) | K (eb / PRE) | PRE coverage |
|---|---|---|---|---|---|
| c_D 10⁻² | 0 / 1 | 0 / 0 | 7.7 / 28.3 / 1.0 | 113 / 2,723 | 48% |
| c_D 3×10⁻³ | 0 / **6** | 0 / 2 | 8.8 / declines all / 1.0 | 92 / 995 | 76% |
| c_D 10⁻³ | 0 / 0 | 0 / 0 | 6.6 / declines all / 2.8 | 77 / 179 | 100% |
| Δp 3×10⁻² | 0 / 0 | 0 / 5 | 3.0 / 3.0 / 1.05 | 444 / 1,447 | 20% |
| Δp 10⁻² | 0 / 0 | 0 / 2 | 7.9 / 7.9 / 3.2 | 224 / 543 | 30% |
| c_L 10⁻¹ | 0 / 0 | 0 / 42 | 1.6 / 1.6 / 1.0 | 2,227 / 5,504 | 90% |
| c_L 3×10⁻² | 0 / 0 | 0 / 49 | 14.2 / 48.7 / 1.0 | 385 / 11,325 | 96% |

- errbudget's misses are out of 608 searches; PRE's are out of 50 designs. eb-3p is errbudget at
  nCellsInCoarsestLevel = 10.
- At three times the budget PRE misses none of 350 designs. Coverage then reaches 86–100% for c_D but stays
  at 38–40% for Δp, and PRE declines 98–100% of c_L designs.

![PRE's band and outcomes](fig_pre_2d1.pdf "**PRE's band rarely contains the exact value for drag and Δp.** (a) Share of PRE fits whose band m ± 2s contains Q*: 608 Arm 1 fits and 50 Arm 2 designs per case. (b) Outcome of the 608 Arm 1 searches per case at z = 2; the numbers are PRE's misses.")

## Why the methods differ

**PRE misses because its band is too narrow.**
- For c_D and Δp, the median |m − Q*| is 3.2–7.7 times s.
- With the white-noise kernel, s measures the scatter of the runs about the fitted polynomial, not the
  extrapolation error. That error is large here because the mesh range is pre-asymptotic: the observed order
  of c_D rises from 0.9 (n_arc 8→12) to 2.5 (48→64).
- On errbudget's runs, PRE's c_D estimate overshoots Q* by a median of +0.12ε, +0.74ε and +2.2ε, which is
  3–6 times the error of the finest converged run. On PRE's own designs it undershoots (−0.05ε to −0.47ε).
- The selector then accepts runs close to m but more than ε from Q*; the worst miss is 1.84ε.

**PRE declines for two reasons.**
- *c_L.* Its band is about right (70–96% coverage) but wide. At residualControl 10⁻², the lift error is 4–79%
  of its value and changes sign with the mesh, so 2s ≤ ε holds in only 55% and 32% of Arm 1 fits.
- *c_D at 3×10⁻³ and 10⁻³.* The band is narrow but off-centre, so few seen runs qualify.

**PRE's settings are cheaper.** errbudget keeps a margin: its settings have true errors of 0.08–0.36ε
(median). PRE takes the cheapest run inside its band, which is cheaper than errbudget's in 72–100% of
returns. This is Table 7's trade-off: a tighter selection is cheaper and less reliable.

**K differs only through its definition.** Both methods spend the same; PRE's K is larger because K divides
that cost by PRE's cheaper setting.

**errbudget has a limit too.** At nCellsInCoarsestLevel = 10, fine runs that errbudget needs do not converge,
so it declines all 164 searches at c_D 3×10⁻³ and 10⁻³ (eb-3p). PRE returns settings there, and 6 of 50 miss
at 3×10⁻³.

## Verdict

**The pre-stated headline criterion holds in both arms.** Used as a selector, PRE misses targets that errbudget
meets: 8.6% of searches at c_D 10⁻² and 2.0% at 10⁻³ in Arm 1, and 12% at c_D 3×10⁻³ in Arm 2. Its 95% band
contains the exact value in 7–30% of fits for drag and Δp. This supports the paper's argument that a posterior
for Q* is not by itself a certificate for a setting.

**The limits:**
- one steady benchmark with three parameters;
- one PRE variant (SPRE with the white-noise kernel; Matérn and GRE untested);
- a selector of our making, because PRE returns no setting;
- misses near the target (worst 1.84ε);
- PRE's settings are cheaper, and a three-times-larger design removed its misses.

I would present it as one paragraph in the comparison with related work, with the Arm 1 table in the
supplement, not as an abstract-level claim.
