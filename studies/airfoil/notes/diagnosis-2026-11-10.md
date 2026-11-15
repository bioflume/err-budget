# Diagnosis of 10 November 2026: separation, non-convergence, and the free stream

Local runs with OpenFOAM v2512 (`opencfd/openfoam-default:2512`, arm64, in Docker under
Colima on the Mac), driven by `driver.py` unchanged through wrappers that run each OpenFOAM
tool in the container. Hamilton's v2606 behaved identically on the registered set-up (check
and diag runs of 10 November). The runs and scripts are in `results/local/` (not committed):
`run_exp.sh` runs a batch; `code-exp/` is the committed code with the free stream set from
the environment; `code-a3/` is the A3 draft with one extra variant. Every run here is a
diagnosis, not a result of the study.

## 1. The registered set-up fails NACA 0012

Coarsest mesh (n_surf 64, dy1 2×10⁻⁵, R = 50), residualControl 10⁻⁶, α = 4°.

| Set-up | NACA 0012 C_l | m 3 %, t 12 % C_l | Converged |
|---|---|---|---|
| Registered: k∞ 4×10⁻⁷ U², ω∞ 133 U/c, sustaining terms | 0.33 | 0.59 | no: residual 10⁻³–10⁻², lift oscillates |
| Same, no sustaining terms (free stream decays from R = 50) | 0.435 | 0.758 | yes, 1,333 / 1,385 it. |
| ω∞ 5 U/c, μt/μ 0.009, sustaining terms | 0.435 | 0.758 | yes, 1,094 / 1,502 it. |
| ω∞ 5 U/c, μt/μ 1, sustaining terms | 0.435 | 0.758 | yes, 821 / 960 it. |
| ω∞ 5 U/c, k∞ 10⁻⁶ U² (TMR SST-sust), sustaining terms | 0.435 | 0.758 | yes, 822 / 1,075 it. |
| Spalart–Allmaras, ν̃∞ = 3ν (check) | – | 0.759 | yes, 2,704 it. |
| Panel method, inviscid | 0.483 | 0.871 | |

The expected NACA 0012 lift at 4° is about 0.44 (TMR's SST results at Re 6×10⁶ give
C_l = 1.09 at 10°, about 0.109 per degree).

**Cause.** The registered values are TMR's far-field values for standard SST, set hundreds of
chords away and meant to decay before they reach the body. Held at the section by the
sustaining terms, ω∞ = 133 U/c starves the boundary layer: on the upper surface of the
centre design the peak ν_t/ν was 50–63 and fell downstream, against 130–270 with
Spalart–Allmaras (about 170 at x = 0.7 from Clauser's outer-layer eddy viscosity), and the
layer separated at 85 % chord (91 % on NACA 0012). TMR's SST page gives, for SST with
sustaining terms, k_amb = 10⁻⁶ U∞² and ω_amb = 5 U∞/L (L the body length), and warns that
the sustaining terms are not effective with the inflow values usual in external
aerodynamics. With any ω∞ of order U/c the answer is the same to 2×10⁻⁴ in C_l.

## 2. Numerics with the TMR SST-sust free stream

Variants (`af/foamcase.py`): **primary_v1** as registered (SIMPLEC, U 0.9, k ω 0.7,
linearUpwind U, k, ω, limited gradients of k, ω); **upKW** = primary_v1 with upwind k, ω;
**robust** = upKW with U 0.7, k ω 0.5 and a limited gradient of U; **upKWr** = upKW with
U 0.7, k ω 0.5; **fallback_v1** as registered (SIMPLE, p 0.3, U 0.7, k ω 0.7, upwind k ω).
residualControl 10⁻⁶ throughout. "stalls": the largest residual makes no new minimum for
3,000 iterations (level in brackets).

| Design, mesh | primary_v1 | upKW | robust | upKWr |
|---|---|---|---|---|
| m 0, t 8; n64 | 680 | | | |
| m 0, t 12; n64 | | | 822 | |
| m 3, t 12; n64 | 1,075 | | | |
| m 6, t 8; n64 | 950 | | | |
| m 6, t 8; n128, R 500 | | 850 | | |
| m 0, t 16; n128, R 500 | | 568 | | |
| m 6, t 16; n64 | stalls (6.6×10⁻⁵) | 1,670 | 2,010 | |
| m 6, t 16; n64, dy1 4×10⁻⁵ | | stalls (7×10⁻⁵) | 1,417 | 2,112 |
| m 6, t 16; n91 | | stalls | 2,596 | 2,437 |
| m 6, t 16; n91, dy1 10⁻⁵, R 500 | | | stalls, L/D 47 | stalls, L/D 26 |
| m 6, t 16; n128, R 500 | | stalls (4×10⁻⁶), L/D 80.1 | diverges (204) | diverges (144) |
| m 6, t 16; n128, dy1 4×10⁻⁵ | | | diverges (1,080) | diverges (264) |
| m 6, t 16; n181, dy1 10⁻⁵, R 5000 | | see §4 | | |
| m 6, t 16; n256, dy1 10⁻⁵, R 5000 | | see §4 | | |
| m 6, t 14; n64, dy1 4×10⁻⁵ | | stalls | | |
| m 6, t 14; n91, dy1 10⁻⁵, R 500 | | 1,933 | | |
| m 4.5, t 16; n64, dy1 4×10⁻⁵ | | 943 | | |
| m 4.5, t 16; n91, dy1 10⁻⁵, R 500 | | 1,804 | | |

Also: limiting the gradient of U alone (on primary_v1) stalls at m 6, t 16; one
non-orthogonal corrector on robust converges (1,890); SIMPLE with p 0.3 (fallback_v1's
coupling) diverges within 21 iterations even with the new free stream. In the diverging
robust run, k grows from 0.04 to 12 within ten iterations in the first wall cells of the
upper surface (x = 0.33, 0.45, 0.74) before the pressure solver fails.

**Reading.** Only the 6 %-camber thick sections (t 14 %, 16 %) are difficult: with the
corrected free stream their upper boundary layer separates over the last 4–5 % of chord at
4°, and the separated region makes the steady iteration mesh- and path-dependent: limit
cycles at residuals 4×10⁻⁶–7×10⁻⁵ on some meshes, divergence or a much more separated
state (L/D 26–47 instead of about 73) with lower relaxation. No variant converged every
mesh tried at m 6 %, t 16 %; upKW never diverged.

## 3. Consequences

- The Hamilton numerics sweep (`study.sweep_numerics`, tarball 8d7f58d) runs with the
  registered free stream: not needed.
- Amendment A3 must change the free stream (TMR SST-sust). The numerics and, possibly, the
  design space depend on a decision about the 6 %-camber thick sections (see the report of
  10 November in the chat).

## 4. Runs still going when this note was written (upKW)

At n_surf 128, dy1 2×10⁻⁵, R = 500 the iteration creeps for designs that converge on the
other meshes too: m 4.5 %, t 16 % at 2.4×10⁻⁶ after 3,066 iterations, m 6 %, t 14 % at
3.4×10⁻⁶ after 3,307 (m 6 %, t 16 % stalled there at 4×10⁻⁶). m 6 %, t 16 % at n_surf 181
and 256 (dy1 10⁻⁵, R = 5000): largest residuals about 9×10⁻⁶ and 3×10⁻⁶ after 1,672 and
3,688 iterations.
