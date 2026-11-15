# Airfoil use case: pre-registered protocol

**Status.** Registered on 10 November 2026, before any run is submitted. This file is
committed on its own; the commit that adds it is the registration. Later changes are
dated amendments in §14, each made before the data it affects exist. The only changes
allowed after data exist are the mechanical rules of §4.2–4.4, which are fixed here.

**What the study is for.** The article on errbudget needs one realistic demonstration: a
user who runs a parameter study with a black-box code tunes the numerical setting once,
at the hardest member of the study, and reuses it for every member. Here the code is
OpenFOAM, the study is a 25-design aerofoil optimisation, and the observable is the
lift-to-drag ratio.

## 1. Questions

- **Q1 (accuracy).** At the hardest design, does the setting that errbudget returns meet
  the target in L/D from every starting setting on the 540-setting grid? How often does
  it decline, and why?
- **Q2 (practice).** On the same stored runs, how do the procedures of current practice
  fare: the sequential convergence study (S), the sequential study with a joint
  refinement (SJ), both at thresholds η = ε, ε/3 and ε/10, and selection against a
  reference run (Ref.)?
- **Q3 (reuse).** In the 25-design optimisation, how close to the reference optimum, in
  design and in L/D, are the optima found with errbudget's setting, with the setting of
  the OpenFOAM tutorial and with a cautious setting? At what total cost?
- **Q4 (payback).** After how many optimisation runs does errbudget's search pay for
  itself, against each of the other settings?

## 2. Flow, observable and cost

| | |
|---|---|
| Code | OpenFOAM v2606 (ESI), `simpleFoam`, through its input files only |
| Flow | steady incompressible RANS, k-ω SST (`kOmegaSST`), Re = 3×10⁶ on chord c = 1, U∞ = 1, ν = 1/(3×10⁶), α = 4° |
| Sections | *Amended (§14, A4): camber m ∈ {0, 1, 2, 3, 4} %, thickness t ∈ {8, 9, 10, 11, 12} %; corners (0, 8), (0, 12), (4, 8), (4, 12) %, centre (2, 10) %.* As registered: NACA 4-digit, maximum camber at 0.4c; camber m ∈ {0, 1.5, 3, 4.5, 6} %, thickness t ∈ {8, 10, 12, 14, 16} %: 25 designs |
| Trailing edge | closed (sharp): last thickness coefficient −0.1036 instead of −0.1015, so the five coefficients sum to zero. The standard section has a blunt edge of 0.021 t (0.25 % c at t = 12 %) |
| Incidence | the section is pitched 4° nose-up about its trailing edge; the free stream is along x and the wake cut along the x-axis |
| Free-stream turbulence | *Amended (§14, A3): k∞ = 10⁻⁶ U∞² (Tu = 0.082 %), ω∞ = 5 U∞/c (μt/μ = 0.6), TMR's values for SST with the sustaining terms.* As registered: NASA Turbulence Modeling Resource values for SST: k∞ = 9×10⁻⁹ a∞², μt/μ = 0.009, at the TMR NACA 0012 Mach number 0.15: k∞ = 4.0×10⁻⁷ U∞² (Tu = 0.052 %), ω∞ = 133.3 U∞/c. The sustaining terms of Spalart & Rumsey (`decayControl yes; kInf; omegaInf;`) keep them from decaying, so the turbulence that reaches the section does not depend on the far-field radius |
| Wall | wall-resolved (low-Re): `noSlip`, ν_t = 0 (`nutLowReWallFunction`), k = 0, ω from `omegaWallFunction` (viscous-sublayer value below y⁺ = 11) |
| Far field | `freestreamVelocity` and `freestreamPressure` on the whole outer boundary; `inletOutlet` for k and ω |
| Numerics (primary) | *Amended (§14, A3): `bounded Gauss upwind` for k and ω; the rest unchanged.* SIMPLEC (`consistent yes`), relaxation U 0.9, k and ω 0.7; `bounded Gauss linearUpwind` for U, k and ω, with `cellLimited Gauss linear 1` gradients for k and ω; `Gauss linear corrected` Laplacians; no non-orthogonal correctors; `meshWave` wall distance. p: GAMG with GaussSeidel, relTol 0.1; U, k, ω: smoothSolver, GaussSeidel, 2 sweeps, relTol 0.1 |
| Observable | L/D = C_l/C_d from OpenFOAM's `forceCoeffs` function object at the last iteration (lRef 1, Aref 1, magUInf 1, lift along y, drag along x). C_l, C_d and C_m (quarter chord) are recorded |
| Cost | cells × iterations (primary: deterministic, independent of the node). Core-seconds (ExecutionTime × cores) and wall time are recorded |
| Converged | `simpleFoam` reports "SIMPLE solution converged": the initial residuals of p, U, k and ω are all below residualControl |

A run that is not converged is one of: stopped at the iteration cap (20,000; 40,000 for the
reference tolerance 10⁻⁷); stopped as **stalled**, when the largest initial residual has
not reached a new minimum (by 5 %) for 3,000 iterations; or **failed**, when the solver
diverged or crashed after its first iteration. Every such run is a failed run for errbudget
and for the three procedures of practice, as in the article (version 1.2.1 treats an
unconverged run as failed), and it is recorded with its cost. A crash before the first
iteration, a mesh that blockMesh rejects, or a missing module is an infrastructure error,
not a result: it is fixed and the run repeated.

## 3. Numerical parameters

| Parameter | Values, coarse → fine | Class | Refinement scalar h | Model |
|---|---|---|---|---|
| `n_surf` (mesh level): cells along each side of the section | 64, 91, 128, 181, 256 | discretization | 1/n_surf | power |
| `dy1`: first-cell height at the wall (chords) | 4×10⁻⁵, 2×10⁻⁵, 1×10⁻⁵ | discretization | dy1 | power |
| `radius`: far-field radius R (chords) | 50, 500, 5000 | discretization | 1/R | power |
| `res_ctrl`: SIMPLE residualControl of p, U, k, ω | 10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶ | nested solver | value | power |
| `gamg_tol`: absolute tolerance of the GAMG p solver; U, k, ω use 1/100 of it | 10⁻⁶, 10⁻⁷, 10⁻⁸ | nested solver | value | power |

5 × 3 × 3 × 4 × 3 = 540 settings. The user's default start is the middle of each list:
`n_surf` 128, `dy1` 2×10⁻⁵, R 500, `res_ctrl` 10⁻⁵, `gamg_tol` 10⁻⁷.

### 3.1 The mesh and why these values

*Amended (§14, A1 and A2): the rays, the split points, the mesh-quality figures and the
outer layer of this section, and the mesh is computed point by point and written as an
OpenFOAM polyMesh, not by blockMesh; the parameters and their values are unchanged.*

The mesh is a structured C-grid of 16 blocks written by `af/cgrid.py` as a
blockMeshDict, from the NACA formulas:

- an **inner layer** 1.5 chords thick around the section and the near wake (to x = 2.5c),
  whose grid lines leave the wall along straight rays: the surface normal ahead of x/c = 0.3,
  then turning in arc length to the bisector of the trailing-edge normal and the wake normal;
- an **outer layer** from there to a circle of radius R about (0.5c, 0) ahead of x = 0.5c,
  the lines y = ±R behind it, and an outlet at x = 0.5c + R.

The inner layer and the near wake do not depend on R, so the far-field parameter changes
only the cells it adds. All side edges of the inner layer end at the same cell size and all
side edges of the outer layer start from it; otherwise the rows of the thin, tall cells
behind the trailing edge drift and their faces approach 90° of non-orthogonality.

- **Mesh level.** With r = n_surf/64, every surface and wake spacing scales with 1/r and
  every growth ratio g with g^(1/r): each step of the list refines by √2 in each direction,
  except the first-cell height and the radius, which are parameters of their own. Level 0
  has leading- and trailing-edge spacings of 2×10⁻³c and growth ratios 1.25 (normal) and
  1.20 (wake). At the middle `dy1` and R the meshes have 18,000, 34,848, 66,690, 129,130
  and 251,430 cells; over the other two discretization parameters the coarsest spans
  12,768–23,290 and the finest 177,684–336,192. The coarsest still converges (Stage 0
  checks) and is a few per cent wrong; the finest is about as fine as the 897×257 grid of
  the NASA TMR NACA 0012 family and costs under an hour per run (§12).
- **First-cell height.** The flat-plate skin friction at Re = 3×10⁶ (c_f = 0.026 Re^(−1/7))
  gives u_τ = 0.039 U∞, so the cell-centre y⁺ of the three values is nominally 2.4, 1.2 and
  0.6, with peaks two to three times higher near the suction peak. The coarsest is what
  many users accept for a wall-resolved run; the finest is a careful choice.
- **Far-field radius.** Without a point-vortex correction, a finite domain imposes a
  downwash of order Γ/(2πR) at the section, which tilts the force and adds a drag of about
  C_l Δα. Our estimate (unmeasured) for the most cambered sections is an L/D error near
  4/R: 8 % at 50 chords, 0.8 % at 500 and 0.08 % at 5000. Because the outer layer is graded
  geometrically, a tenfold radius costs a few dozen more rows, not ten times the cells. The
  list runs from the 20–50 chords of common practice, through the 500 of the TMR grids, to
  5000.
- **residualControl.** From a loose 10⁻³ to 10⁻⁶; it contains the tutorial's 10⁻⁵.
- **GAMG tolerance.** From the tutorial's 10⁻⁶ to 10⁻⁸. No value exceeds the smallest
  residualControl: a linear tolerance above residualControl stops solving p before the outer
  iteration can converge, so those runs would all reach the cap, at the price of the cap.
  The smoothSolver tolerance of U, k and ω is 1/100 of it, the tutorial's ratio.

The geometric check of `af/cgrid.py` (transfinite interpolation of the block edges, as
blockMesh builds the interior) finds no inverted cell, at most 25° of non-orthogonality
and at most 13° between a wall ray and the wall normal, on every design at every level.

### 3.2 Far-field check after Stage 0 (mechanical)

*Withdrawn (§14, A5): the radius list stays 50, 500, 5000; the Stage-0 analysis reports the
far-field error at R = 5000 instead.*

If at any corner design of Stage 0 (runs D, §4.1) the L/D at R = 5000 differs from that at
R = 50,000 by more than 0.1 % (ε/10), the radius list becomes 500, 5000, 50,000 and the
reference radius 500,000, for Stages 1 and 2. The Stage-0 analysis writes the factor into
`config/stage1.json`.

### 3.3 Numerics check after Stage 0 (mechanical)

*Amended (§14, A3): the fallback numerics are withdrawn. If the rule below triggers, the
Stage-0 analysis prints STOP, and the study stops for an amendment instead of switching.*

If any Stage-0 run diverges, or more than 2 of its 17 runs at residualControl 10⁻⁶ stop
unconverged, every later run uses the **fallback numerics**: SIMPLE (`consistent no`) with
the airFoil2D relaxation (p 0.3, U 0.7, k and ω 0.7) and first-order `bounded Gauss upwind`
for k and ω. Stage 0 is then repeated with the fallback, and the repeat replaces it for
§4.2 and the cost model; both are reported. Fallback runs carry the tag suffix `_fb`.

## 4. Stages

### 4.1 Stage 0 (30 runs, partition `shared`)

*Amended (§14, A4): the centre design is m 2 %, t 10 % and the corners those of the A4 design
space; the runs are otherwise as below.*

Middle setting unless stated.

- **A.** Centre design (m 3 %, t 12 %) at n_surf 64, 128, 256: time per run against cells,
  iterations to convergence, y⁺ achieved, forces written.
- **A′.** Centre design at n_surf 256 with the finest tolerances (10⁻⁶, 10⁻⁸): the slowest
  convergence on the Stage-1 grid.
- **B.** Centre design at n_surf 181, on one core and on four: parallel speed-up.
- **C.** The four corner designs (m, t) ∈ {0, 6} % × {8, 16} % at n_surf 128 and 181:
  convergence, and the hardest design (§4.2).
- **D.** The four corners at n_surf 91, finest tolerances, R = 50, 500, 5000, 50,000: the
  far-field check (§3.2) and an early look at the far-field error. *Amended (§14, A5): R = 50,
  500 and 5000 only, settled; 26 runs in Stage 0.*

A 15-minute check job on the `test` partition runs two tiny cases end to end first.

### 4.2 The hardest design (mechanical)

For each corner, δ = |L/D(n_surf 181) − L/D(n_surf 128)| / |L/D(n_surf 181)|, from runs C.
The hardest design is the corner with the largest δ. A corner whose run at 128 or 181
failed or did not converge ranks above every corner whose runs converged. Corners within
10 % of the largest δ are tied; a tie goes to the thickest, then the most cambered.

### 4.3 Stage 1 (540 + 5 runs)

The hardest design at all 540 settings, plus the reference runs of §5 at that design.
Every run is stored, so that errbudget and the procedures of practice can be searched from
every starting setting on stored runs, as in the article.

### 4.4 Stage 2 (25 designs × 4 settings)

The settings (§7) are written to `config/stage2_settings.json` by the Stage-1 analysis.
Runs already done in Stage 1 (the hardest design at any of the four settings) are reused.

## 5. The reference and its uncertainty

### 5.1 The reference setting (d)

*Amended (§14, A4): a run on a mesh finer than the Stage-1 lists starts from the solution of
the same design at n_surf 91, dy1 2×10⁻⁵, the same radius, mapped onto its mesh; that run is
part of it and its cost is added.*

*Amended (§14, A5): R = 5000 instead of 50,000; the run is settled; the reference value is
L/D(d) plus a far-field correction measured on the n_surf 91 mesh.*

`n_surf` 362, `dy1` 5×10⁻⁶, R 50,000, `res_ctrl` 10⁻⁷, `gamg_tol` 10⁻⁹: one step beyond the
finest value of every list (846,720 cells).

### 5.2 Its uncertainty, and when a target is graded

*Amended (§14, A5): the radius sibling is R = 500; u_adm adds 20 % of the far-field
correction.*

For each design, u = |L/D(d) − L/D(c)| / |L/D(d)|, its change from the cautious setting
(c), as in the article. A target ε is graded at the hardest design only if ε ≥ 10 u there.

Two further runs at the hardest design test this estimate, and it is replaced if they show
that it understates the error of (d):

- three **siblings** of (d), each with one discretization parameter back at its finest list
  value (`n_surf` 256; `dy1` 10⁻⁵; R 5000), show which parameter carries u;
- (d) at **`n_surf` 512** gives the observed order p of the mesh from n_surf 256, 362 and 512
  (the other parameters at their reference values). If p < 2, where the change from (c) can
  understate the mesh error of (d), the admission uses
  u_adm = max(u, |L/D(512) − L/D(362)| / ((√2)^p − 1) / |L/D(d)|).

### 5.3 Targets

ε = 1 % of L/D (primary) and 0.3 % (secondary), relative to the reference value at the
hardest design: ε_abs = ε |L/D(d)|. If 1 % is not admitted, a fallback target of 3 % is
graded and reported as such.

## 6. Grading and outcomes (Stage 1)

The tensor is the 540 stored runs at the hardest design: L/D, cost (cells × iterations)
and the converged flag.

- **The program.** errbudget version 1.2.1, the engine of the article (`errbudget.V121`),
  unchanged, with the classes and refinement scalars of §3. It is replayed from every
  setting whose own run converged.
- **Practice.** S and SJ at η = ε, ε/3, ε/10, taking the parameters in the order of §3
  (discretization first), and Ref.: the cheapest setting within ε of the finest setting of
  the grid. From the same starts. S, SJ and Ref. cannot decline.
- **Outcomes per start**: reached (returned setting converged and within ε of the
  reference); missed (returned and outside ε, or returned an unconverged run); declined
  (program only; split by whether some converged setting on the grid meets ε, and by the
  reason the program gives).
- **Measures**, as in Section 6.3–6.4 of the article: the cost ratio of each returned
  setting to the cheapest converged setting that meets ε; the distinct runs used; the
  search cost K, the summed cost of the runs a search used divided by the cost of one run at
  the setting it returned.
- **The user's case**: the program from the default start, its setting, its error and K.

The replay code is copied from the article's `paper/replay.py`, `baselines.py`, `pooled.py`
and `exp72_baselines_pooled.py`, changed only where it loads data.

## 7. The four settings of Stage 2

*Amended (§14, A6): Stage 2 is not run; the study ends with Stage 1.*

- **(a) errbudget.** The setting errbudget returns from the default start at ε = 1 %, run as
  a user would (target relative to L/D at the start, `relative: true`), found by replaying
  the stored Stage-1 runs, which gives the setting a live run would. If it declines, (a) is
  not run and the decline is the result (§11).
- **(b) tutorial.** The tolerances of OpenFOAM's airFoil2D tutorial
  (`tutorials/incompressible/simpleFoam/airFoil2D/system/fvSolution`, v2606: p GAMG tolerance
  10⁻⁶, relTol 0.1; U and ν̃ smoothSolver 10⁻⁸, relTol 0.1; residualControl 10⁻⁵ for p, U,
  ν̃), that is `res_ctrl` 10⁻⁵ and `gamg_tol` 10⁻⁶, on the middle mesh: `n_surf` 128, `dy1`
  2×10⁻⁵, R 500.
- **(c) cautious.** The finest value of every list.
- **(d) reference.** §5.1.

## 8. The optimisation

For each setting, L/D at the 25 designs is fitted by a full quadratic in (m, t), scaled to
[−1, 1]², by least squares (6 coefficients, 25 points). Its optimum is the maximum of the
fitted surface over the closed box m ∈ [0, 6] %, t ∈ [8, 16] %. The best of the 25 designs is
reported too. The reference optimum is the same construction with (d).

For each of (a)–(c), against (d):

1. the distance of its optimum from the reference optimum, |Δm| and |Δt|;
2. its predicted optimal L/D against the reference optimal L/D;
3. the **loss**: the reference surface at the reference optimum minus the reference
   surface at this setting's optimum, relative: what choosing the design with this setting
   costs in true L/D;
4. whether its best grid design is the reference's;
5. the largest error over the 25 designs of its L/D against (d), relative;
6. the quality of the fit (RMS residual, R²).

## 9. Cost and payback

- The cost of a setting is the sum over the 25 designs of cells × iterations, also in
  core-hours.
- errbudget's search cost is the summed cost of the runs it used from the default start, in
  the same units, and K is that divided by the cost of one run of (a) at the hardest design.
- The total cost of (a) is the search plus its 25 runs; of (b)–(c), their 25 runs.
- With μ the ratio of the mean cost per design of another setting to that of (a), the search
  pays for itself after n* = K/(μ − 1) runs (μ > 1; otherwise never). n* is reported against
  (b), (c) and (d), with the number of optimisation runs (25) beside it.

## 10. What counts as a result worth reporting

The use case is reported in the article whatever it shows; nothing in this protocol depends
on the outcome. The outcome that would support errbudget's claims has all of:

- at an admitted target, no accepted setting misses, from any start;
- setting (a) puts the optimum within one grid step of the reference optimum (|Δm| ≤ 1.5 %,
  |Δt| ≤ 2 %) with a loss below ε;
- the search pays for itself within the 25 runs against (c) (n* ≤ 25).

It is a stronger result if, in addition, the tutorial setting (b) misses the target at the
hardest design, or moves the optimum by more than one grid step, or has a loss above ε.

## 11. Outcomes reported as negative

- Any accepted setting that misses an admitted target, with its cause.
- errbudget declining from the default start at 1 % (the use case then shows the decline and
  the parameter it names).
- The 1 % target not admitted (u_adm > 0.1 %): reported as "not graded", with the parameter
  that carries u (§5.2); the 3 % fallback is graded instead.
- The optimum of (a) farther from the reference optimum than that of (b) or (c), or a loss
  of (a) above ε.
- n* > 25 against (c): the search does not pay for itself in this study.
- A stage stopped by its cost guard (§12).

## 12. Budget and stop rules

Predicted at 3 µs per cell and iteration (COST.md; unmeasured until Stage 0): Stage 0, 3
core-hours; Stage 1, 60 plus 20 for the reference runs; Stage 2, about 120; about 200 in
all. The agreed budget is 600 core-hours and the ceiling 1,000.

- Guards: Stage 0, 8 core-hours; Stage 1 with its reference runs, 400; Stage 2, 500. A stage
  that passes its guard is stopped, and what remains is reported.
- If the cost model measured in Stage 0 predicts more than 900 core-hours for Stages 1 and 2,
  the study stops before Stage 1 and the design is revisited (an amendment).

## 13. Software

- errbudget engine 1.2.1 from the article's repository (commit 4d32148), unchanged.
- This study: branch `study/airfoil`, directory `airfoil/`. Each run's JSON records the hash
  of the driver files, the OpenFOAM version and the Slurm job.

## 14. Amendments

**A1 (10 November 2026, before any valid run).** On Hamilton the two runs of the check job
diverged after 46 and 83 iterations, and checkMesh reported a maximum non-orthogonality of
88.7° on the coarsest mesh, where the generator's own check had reported at most 25°. The
check had rebuilt each block's interior by transfinite interpolation and measured cell
centres as corner means; blockMesh interpolates differently (`block::createPoints`), and
OpenFOAM uses centroids. On the thin, tall cells that the trailing-edge clustering carries
out to the far field (aspect ratios of 10⁴–10⁶), a taper of a fraction of a per cent along a
column moves the centroid enough to give faces near 90°. Changed, in the mesh family of §3.1
only (not the parameters, their values, the targets or any rule):

1. The generator's check reproduces blockMesh's interpolation and OpenFOAM's centroids. On
   the first generator's two check meshes it gives checkMesh's 88.70° and 88.46° and its
   minimum cell volumes (a test keeps it so).
2. The split points are where the surfaces cross x = 0.5c, the centre of the far-field
   circle, and their rays are vertical. Behind them the mid boundary is the line y = ±1.5c,
   where the rays end, so the outer rear and wake blocks are exact rectangles.
3. Ahead of the split points the rays start as the leading-edge normal and turn smoothly to
   the vertical, never following a concave stretch of surface; behind them they turn to the
   trailing-edge bisector and stay parallel over the last 0.15c, and in the near wake they
   stay parallel for 0.15c before turning to the vertical.

The rays now differ in length, and blockMesh's interpolation tilts them a little near the
wall: the first-cell height normal to the wall is 0.79–1.06 dy1 along the surface (exactly
dy1 on the block edges), and scales with dy1. The cell counts are unchanged. Over all 25
designs at seven mesh specifications from the coarsest to the reference, and the four corners
at all 75 combinations of level, first cell and radius (the amended radii of §3.2 included):
no inverted cell, no face above 70°, at most 52° (two cells behind the trailing edge on the
coarsest mesh with the finest wall spacing) and at most 17° at the wall. The two check runs,
made with the first mesh, are discarded.

**A2 (10 November 2026, before any valid run).** The mesh is computed point by point in
Python and written directly as an OpenFOAM polyMesh; blockMesh is no longer used.

*Reason.* In the blockMesh grid of A1 the columns over the rear of the section kept the
width of the surface cells all the way to the far field, while the columns ahead of the
split points fanned out with the radius. Where the two met, at x = 0.5c next to the far
boundary, neighbouring cells differed in size by a factor that grows with R: 36 at R = 50,
364 at R = 500, 3,543 at R = 5000, and about 35,000 at the reference radius. blockMesh fills
a block by straight-line interpolation between its edges and cannot give the curved,
widening columns that remove the jump.

*Evidence.* The diagnostic sweep of 10 November (five runs at R = 500 or on the coarsest
mesh, both numerics variants of §3.3, the 6 % camber corners): every run at R = 500 diverged
within 17–20 iterations, with either numerics; the check runs at R = 50, where the jump was
36, ran (one converged, the m 6 %, t 16 % corner oscillated). These runs were made with the
A1 mesh and are discarded with it.

*The mesh now* (the parameters, their values and the cell counts are unchanged):

1. The inner layer is placed explicitly along the rays of A1, so the first cell is exactly
   dy1 along every ray (0.97–1.00 dy1 normal to the wall). On the wake cut the first cell
   grows from dy1 at the trailing edge to the wake spacing smoothly over the near wake.
   The front rays end on the convex hull of their ends and the lines y = ±1.5c, so the mid
   boundary is convex.
2. The outer layer's rows are offset curves of the mid boundary, at distances growing
   geometrically to R − 1.5c: the far field is the line y = ±R behind x = 0.5c and, ahead,
   a curve about R from the leading edge; the outlet is at x = 0.5c + R. Its columns start
   along the normals of the mid boundary and lean, in proportion to the distance, to a far
   distribution spread evenly along the far boundary (85 % even, 15 % along the normals).
3. Over all 25 designs at seven mesh specifications from the coarsest to the reference, and
   the four corners at all 75 combinations of level, first cell and radius (the radii of
   §3.2 included) and at n_surf 512: no inverted cell, no face above 70°, at most 57°
   (a few faces where the outer layer starts behind the trailing edge, on the coarsest
   meshes; below 50° from n_surf 128), neighbouring cells differing in size by at most a
   factor of 6.4, and at most 13° at the wall.

*Checks.* checkMesh runs on every mesh. A run whose mesh fails a check is an error, not a
result, with one exception: checkMesh's high-aspect-ratio check, which fails above a fixed
ratio of 1000 that every wall-resolved mesh at y⁺ ≈ 1 exceeds (4,673 on the coarsest mesh).
Severely non-orthogonal faces (above 70°), which checkMesh reports as a warning, also fail
a run. Before the diagnostic or Stage-0 runs, the check job runs checkMesh on the coarsest
and finest meshes of the Stage-1 lists for the corners m 6 %, t 16 % and m 0, t 8 %, on the
reference mesh and on the n_surf 512 mesh, and fails if any fails. The writer is
unit-tested: the face count, owner below neighbour in OpenFOAM's upper-triangular order,
outward normals on the boundary, closed cells (the face area vectors of each cell sum to
zero), positive volumes, and the files read back. Results carry a mesh version; results of
an older mesh are set aside, not used.


**A3 (10 November 2026, before any valid run).** The free-stream turbulence of §2, the
convection of k and ω (§2, numerics), and the fallback of §3.3.

*Reason.* With the A2 mesh, every run with the registered set-up behaved alike, on Hamilton
(v2606) and in a local OpenFOAM v2512 (Docker): the primary numerics never converged (the
largest initial residual stayed at 10⁻³–10⁻², lift oscillating), the fallback diverged
within 50 iterations, and at 4° the upper boundary layer separated at 75–91 % chord, with
lift 25–40 % below the inviscid value of a panel method. NACA 0012 (m 0, t 12 %, a design
of the study) gave C_l = 0.33; about 0.44 is expected (TMR's SST results at Re 6×10⁶,
C_l = 1.09 at 10°). The cause is the free-stream turbulence of §2. TMR's values for
standard SST (ω∞ = 133 U∞/c) are set at a far field hundreds of chords away and decay
before they reach the body. Held at the section by the sustaining terms, they starve the
boundary layer: on the upper surface of the centre design the peak ν_t/ν was 50–63 and fell
downstream, against 130–270 with the Spalart–Allmaras model on the same mesh, and the layer
separated. TMR's page for SST with the sustaining terms gives k_amb = 10⁻⁶ U∞² and
ω_amb = 5 U∞/L, with L the body length, and warns that the sustaining terms are not
effective with the inflow values usual in external aerodynamics.

*Evidence* (local runs, `notes/diagnosis-2026-11-10.md`, commit 5b2cd02; coarsest mesh,
residualControl 10⁻⁶). With the registered free stream: NACA 0012 C_l 0.33 and the centre
design 0.59, neither converged. With k∞ = 10⁻⁶ U∞², ω∞ = 5 U∞/c: 0.435 and 0.758, both
converged (822 and 1,075 iterations), attached; Spalart–Allmaras gives 0.759 for the centre
design; other free streams with ω∞ of order U∞/c, or none sustained, give the same C_l to
2×10⁻⁴. With the new free stream the registered primary numerics converge at m 0, t 8;
m 3, t 12; and m 6, t 8 (680–1,075 iterations), but at m 6 %, t 16 % the largest residual
stops at 6.6×10⁻⁵: that section separates over its last 4–5 % of chord at 4°. First-order
convection of k and ω converges there (1,670 iterations), at m 0, t 16 % and m 6 %, t 8 %
on n_surf 128 with R = 500 (568 and 850 iterations), and at m 4.5 %, t 16 % on two coarser
meshes; it never diverged. But on several meshes of the thick cambered sections it stalls:
at m 6 %, t 16 % on the coarsest mesh with dy1 4×10⁻⁵ and on n_surf 91 (largest residual
6–7×10⁻⁵), and at m 6 %, t 16 %, m 6 %, t 14 % and m 4.5 %, t 16 % on n_surf 128 with
R = 500, where the largest residual levels off at 2–4×10⁻⁶; at m 6 %, t 16 % on n_surf 181
and 256 with R = 5000 it stood at 3–8×10⁻⁶ after about 2,000 and 4,000 iterations. Variants with lower relaxation (U 0.7, k and ω 0.5, with
or without a limited gradient of U) converged some of those meshes but diverged on others
(k running away in the first wall cells) or settled in a much more separated state. SIMPLE
with p relaxed at 0.3, the coupling of the registered fallback, diverged within 21
iterations even with the new free stream.

*Changed* (the parameters, their values, the targets, the designs, the mesh and every other
rule are unchanged):

1. Free-stream turbulence (§2): k∞ = 10⁻⁶ U∞² (Tu = 0.082 %), ω∞ = 5 U∞/c (μt/μ = 0.6),
   with the sustaining terms.
2. Primary numerics (§2): `bounded Gauss upwind` for k and ω instead of `bounded Gauss
   linearUpwind`; everything else as registered. The convection of U stays second order.
3. Fallback (§3.3): withdrawn. If the rule of §3.3 triggers, the study stops and the cause
   is reported; a further amendment decides how to go on.
4. A2's figure for the aspect ratio of the coarsest mesh is checkMesh's own, 1,976 (4,673
   was the generator's ratio of edge lengths); the exception of A2 stands.
5. Results carry a set-up version; results made with the registered free stream or numerics
   (the check, diagnostic and numerics runs of 10 November) are set aside, not used. The
   study's code remains OpenFOAM v2606 on Hamilton: v2512 served only this diagnosis, and
   the check job is repeated on Hamilton before any further run.

*Before Stage 0* (diagnostic, not a stage): a convergence map at the 6 %- and 4.5 %-camber
thick sections (t 14, 16 %): the 45 Stage-1 meshes at residualControl 10⁻⁶, GAMG 10⁻⁸, and
the reference setting at m 6 %, t 16 % and m 4.5 %, t 16 %, with the numerics of this
amendment. Whether the design space, or the tolerance lists, stay as registered is decided
on it (in particular, whether the reference reaches residualControl 10⁻⁷), by a further
amendment if they change. Its runs are stored like any other; those that coincide with a
later run of a stage are reused by it.


**A4 (11 November 2026, before any valid run).** The design space (§2), and the start of the
runs on meshes finer than the Stage-1 lists (§5.1).

*Reason.* The convergence map of A3 (`notes/convmap-2026-11-11.md`, commit 708dc0a; 182 runs on
Hamilton) found the thick cambered sections unusable at 4°: at residualControl 10⁻⁶ the 45
Stage-1 meshes converged 23 (m 6 %, t 16 %), 28 (6, 14), 32 (4.5, 16) and 39 (4.5, 14) times,
and in the stalled runs L/D still moved by about 1 % over their last 1,000 iterations (the
targets are 1 % and 0.3 %). Both reference runs diverged within 150 iterations.

*Evidence* (local runs, OpenFOAM v2512; `notes/` and the chat record of 11 November):

1. The stalls come from the flow, not the mesh or the turbulence model. On the meshes where
   all four thick designs stalled, m 3 %, t 12 % and NACA 0012 converged (380–1,920
   iterations). In the stalled runs the solution keeps changing in the first wake cells
   behind the trailing edge, where the upper boundary layer has separated (over 0.2 % of the
   upper surface at m 4.5, t 14 % to 5 % at m 6, t 16 %). The Spalart–Allmaras model stalled
   on the same meshes (6,900–19,000 iterations), and the airFoil2D tutorial's SIMPLE
   coupling diverged with it.
2. On the four meshes where the thick designs failed most, m 4 %, t 12 % converged on all
   four (845–2,230 iterations); m 6 %, t 12 % and m 6 %, t 16 % at 2° on two of four.
3. The reference diverges from the uniform start: m 4 %, t 12 % at R = 50,000 after 172
   iterations, as m 6 %, t 16 % and m 4.5 %, t 16 % did on Hamilton. A potential-flow start
   (13–24 iterations), lower relaxation (46) and first-order convection of U (200) failed
   too. A far-wake mesh that caps the cells' aspect ratio was tried and abandoned: with it
   the reference diverged after 107 iterations, and the finest Stage-1 mesh (n_surf 256,
   dy1 10⁻⁵, R = 5,000) diverged after 109, which on the A2 mesh does not diverge. Started
   from the converged n_surf 91 solution of the same design and radius, mapped onto the
   reference mesh, the reference of m 4 %, t 12 % converged to residualControl 10⁻⁷ in
   10,316 iterations (L/D 75.890; it still moved by 0.16 % over the last 1,000); the
   n_surf 91 run cost 3 % of it.

*Changed* (the mesh of A2, the free stream and numerics of A3, the parameters and their
lists, the reference setting, the targets and every rule are unchanged):

1. Design space (§2): camber m ∈ {0, 1, 2, 3, 4} %, thickness t ∈ {8, 9, 10, 11, 12} % (25
   designs; maximum camber at 0.4c, α = 4° as before). Corners (0, 8), (0, 12), (4, 8),
   (4, 12) %; centre (2, 10) %. The check job runs the centre and the corner m 4 %, t 12 %;
   its mesh check uses the corners m 4 %, t 12 % and m 0, t 8 %.
2. Start of the finest runs (§5.1): a run with n_surf above the Stage-1 list (the reference
   setting, and the n_surf 512 mesh-order check) first runs the same design at n_surf 91,
   dy1 2×10⁻⁵, its own radius, residualControl 10⁻⁶ and GAMG 10⁻⁸ from the free stream, and
   starts from that solution mapped onto its mesh (OpenFOAM's mapFields, internal fields).
   The first run's iterations, cells × iterations and core-seconds are recorded with the
   run and added to its cost. A first run that diverges makes the run an error; one that
   stalls is used as it stands, and recorded. Every other run starts from the free stream.
3. Cost model: 4.6 µs per cell and iteration (median of the convergence map, against 3 µs
   assumed); `COST.md` is recomputed.

*Before Stage 0* (diagnostic, not a stage): a confirmation map on Hamilton: the corners
m 4 %, t 12 % and m 4 %, t 8 % on the 45 Stage-1 meshes at residualControl 10⁻⁶ and GAMG
10⁻⁸, and the reference of m 4 %, t 12 % started as in item 2. Locally, the finest Stage-1 mesh of
m 4 %, t 12 % (n_surf 256, dy1 10⁻⁵, R = 5,000) stalled at residualControl 10⁻⁶ after 5,018
iterations (L/D 81.38). If many runs do not converge at the tightest tolerances, the
tolerance lists go back to the user, and any change is a further amendment.


**A5 (13 November 2026, before any valid run).** The reference (§5.1, §5.2), Stage 0's
far-field runs (§4.1, D) and the far-field check (§3.2).

*Reason.* The reference of the confirmation map of A4 (`notes/confirm-2026-11-12.md`) stalled
after 10,978 iterations with L/D in a slow, damped oscillation of about ±2 % while its
residuals approached 10⁻⁷. On the reference mesh of m 4 %, t 12 %, run for 40,000
iterations without a residual stop (`notes/settle-2026-11-13.md`), L/D at R = 50,000 never
settled to 10⁻⁴ over 2,000 iterations (to 10⁻³ only after 34,294; 77.64 at 10,000 and 77.96
at 40,000 iterations) although its residuals fell below 10⁻⁷; at R = 5000 it settled from
iteration 8,626 (L/D 80.368). The residual criterion leaves L/D moving more the larger the
domain: in the confirmation map, at residualControl 10⁻⁶, the median spread of L/D over the
last tenth of the iterations was 0.02 % at R = 50, 0.1 % at 500 and 0.2 % at 5000. A reference
at R = 50,000 cannot be converged to within the targets at an affordable cost, and the check
of §3.2, which compares R = 5000 with R = 50,000 to 0.1 %, cannot work. The far-field error at
R = 5000 is small: fits L/D = Q + A/R to R = 50 and 500 give A between −350 and −430 on every
mesh of m 4 %, t 12 %, about 0.1 % of L/D at R = 5000.

*Changed* (the parameters and their lists, the mesh, the set-up, the targets and every other
rule are unchanged):

1. **Settled runs.** A run that must be accurate stops only when its largest initial residual
   is at its residualControl and the relative spread of L/D over the preceding 2,000
   iterations is at most 10⁻⁴. OpenFOAM's own residual stop is switched off for it
   (residualControl 10⁻¹² in its dictionary); the driver applies both conditions; the
   iteration cap stays (40,000 at residualControl 10⁻⁷, otherwise 20,000), and a run that
   reaches it unsettled is unconverged. Settled runs carry the tag suffix `_st`. Settled are:
   the reference setting (d), its siblings and the n_surf 512 run of §5.2; the first runs of
   A4, item 2; the Stage-0 runs D; and the runs of item 3.
2. **The reference setting (§5.1):** n_surf 362, dy1 5×10⁻⁶, residualControl 10⁻⁷, GAMG
   10⁻⁹ as before, at R = 5000, the finest listed radius, instead of 50,000; settled; started
   as in A4.
3. **Far-field correction.** For every design with a reference run, the first run of A4
   (n_surf 91, dy1 2×10⁻⁵, R = 5000, residualControl 10⁻⁶, GAMG 10⁻⁸, settled) and one more
   run, the same at R = 500 (settled), give Δ_R = (L/D(5000) − L/D(500)) / 9, the first-order
   extrapolation in 1/R from R = 5000 to an unbounded domain. The reference value of the
   design is L/D(d) + Δ_R wherever the protocol uses L/D(d): u, the targets, the grading and
   Stage 2.
4. **Its uncertainty (§5.2).** u_adm = max(u, u_512) + 0.2 |Δ_R| / |L/D(d)|, where u_512 is
   the n_surf 512 term of §5.2: the correction is taken as uncertain by 20 %, twice the spread
   of A across the mesh levels of the confirmation map. The radius sibling of (d) is R = 500,
   the next value down.
5. **Stage 0, runs D (§4.1):** the four corners at n_surf 91, finest tolerances, at R = 50, 500
   and 5000 only, settled: 12 runs, 26 in Stage 0; the numerics check of §3.3 counts the 13
   Stage-0 runs at residualControl 10⁻⁶.
6. **§3.2 is withdrawn.** The radius list stays 50, 500, 5000 and the reference radius 5000.
   The Stage-0 analysis reports, for each corner, A and Δ_R from the runs D.
7. **Cost.** A settled reference run takes about 11,000 iterations (8,626 to settle on the
   reference mesh at R = 5000) and the correction run about 1.6 core-hours per design.
   Stages 0–2 are predicted at about 405 core-hours central (`COST.md`).

**A6 (14 November 2026, after Stage 1).** Stage 2 (§7–§9) and the end of the study. Decided
after the Stage-1 results, with the user's approval.

*Reason.* Stage 1 (`notes/stage1-2026-11-14.md`): the reference stalled with its L/D settled
but its pressure residual just above 10⁻⁷, and the n_surf 512 run diverged. With the
reference's L/D taken as information only, u = 4.5 % and no target is admitted, the 3 %
fallback included; the reference's dy1 sibling differs from it by 2 %. errbudget declined from
the default start at 1 %, so (a) is not run (§7). Stage 2 would compare only (b) and (c) with a
reference uncertain by about 2 %, and would answer neither Q3 nor Q4 for errbudget.

*Changed:*

1. **Stage 2 is not run.** Q3 and Q4 are not answered by this study and are reported as such
   (§11).
2. **Stage 1 stands as run.** Neither the reference nor the n_surf 512 run is rerun. Reported:
   the accounting, the outcome of the reference runs, errbudget's decline from the default
   start and the parameter it names, and, as information only, the grading of §6 against the
   stalled reference's L/D.
3. **The study ends here.** No further Hamilton runs are made for it.
