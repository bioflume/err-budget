"""The study as data: flow, designs, parameters, settings and sweeps (PROTOCOL.md).

Everything that PROTOCOL.md fixes and the code needs lives here, so the driver,
the job file, the Slurm scripts and the analysis read one definition.
"""
from __future__ import annotations

import itertools
import json
import math
import os
from dataclasses import dataclass, asdict

import numpy as np

# ---------------------------------------------------------------- the flow
RE = 3.0e6
U_INF = 1.0
CHORD = 1.0
NU = U_INF * CHORD / RE
ALPHA_DEG = 4.0
P_CAMBER = 0.4
# free-stream turbulence (PROTOCOL.md 14, A3): the NASA Turbulence Modeling Resource's
# values for SST with the sustaining terms (SST-sust), k_amb = 1e-6 U^2 (Tu = 0.082 %),
# omega_amb = 5 U / c (mu_t/mu = 0.6).  The values first registered, TMR's far-field
# values for standard SST (k 4.0e-7 U^2, omega 133 U/c), are meant to decay over hundreds
# of chords; held at the section by the sustaining terms they starved the boundary layer.
K_INF = 1e-6 * U_INF**2
OMEGA_INF = 5.0 * U_INF / CHORD
# version of the case set-up apart from the mesh (free stream, numerics): results made
# with an older set-up are set aside like results of an older mesh.  1: as registered;
# 2: amendment A3
SETUP_VERSION = 2

# ---------------------------------------------------------------- designs
# amendment A4 (PROTOCOL.md 14): camber 0-4 %, thickness 8-12 %.  As registered: camber
# 0, 1.5, 3, 4.5, 6 % and thickness 8-16 %, whose thick cambered sections separate at the
# trailing edge at 4 degrees and do not converge on many meshes (convergence map, 11 Nov)
CAMBERS = [0.0, 0.01, 0.02, 0.03, 0.04]
THICKNESSES = [0.08, 0.09, 0.10, 0.11, 0.12]
DESIGNS = [(m, t) for m in CAMBERS for t in THICKNESSES]
CORNERS = [(0.0, 0.08), (0.0, 0.12), (0.04, 0.08), (0.04, 0.12)]
CENTRE = (0.02, 0.10)

# ---------------------------------------------------------------- parameters
# values run coarse (cheap) -> fine; errbudget's default start is the middle of each
PARAMS = [
    # name        values                              kind              h
    ("n_surf",   [64, 91, 128, 181, 256],             "discretisation", "1/value"),
    ("dy1",      [4e-5, 2e-5, 1e-5],                  "discretisation", "value"),
    ("radius",   [50.0, 500.0, 5000.0],               "discretisation", "1/value"),
    ("res_ctrl", [1e-3, 1e-4, 1e-5, 1e-6],            "nested-solver",  "value"),
    ("gamg_tol", [1e-6, 1e-7, 1e-8],                  "nested-solver",  "value"),
]
NAMES = [p[0] for p in PARAMS]
LADDERS = [p[1] for p in PARAMS]
KINDS = [p[2] for p in PARAMS]
START_IDX = tuple(len(v) // 2 for v in LADDERS)          # (2, 1, 1, 2, 1)
GRID = list(itertools.product(*[range(len(v)) for v in LADDERS]))   # 540 settings


def hs_of(ladders):
    """Refinement scalars of the given lists, decreasing coarse -> fine (errbudget's h)."""
    out = []
    for (name, _, kind, h), vals in zip(PARAMS, ladders):
        out.append([1.0 / v for v in vals] if h == "1/value" else [float(v) for v in vals])
    return out


def hs():
    return hs_of(LADDERS)


def setting_of_ladders(idx, ladders):
    return {n: ladders[k][j] for k, (n, j) in enumerate(zip(NAMES, idx))}


# The step beyond the finest value of every list (the reference) and the solver caps
# A5: the reference at R = 5000, the finest listed radius (registered: 50,000, where L/D did
# not settle in 40,000 iterations); its value carries a far-field correction (far_field_run)
REFERENCE = dict(n_surf=362, dy1=5e-6, radius=5000.0, res_ctrl=1e-7, gamg_tol=1e-9)
L6 = dict(REFERENCE, n_surf=512)            # mesh-order check of the reference (hardest design)
# OpenFOAM v2606 tutorials/incompressible/simpleFoam/airFoil2D/system/fvSolution:
# p GAMG tolerance 1e-06 relTol 0.1; (U|nuTilda) smoothSolver tolerance 1e-08 relTol 0.1;
# SIMPLE residualControl p 1e-5, U 1e-5, nuTilda 1e-5.  On the middle mesh.
TUTORIAL_TOLERANCES = dict(res_ctrl=1e-5, gamg_tol=1e-6)
DEFAULT_START = {n: v[i] for n, v, i in zip(NAMES, LADDERS, START_IDX)}
MID = dict(DEFAULT_START)
SOLVER_RATIO = 0.01          # smoothSolver (U, k, omega) tolerance = gamg_tol * 0.01, the tutorial's ratio
MAX_ITER = 20000             # iteration cap of every run except the reference
MAX_ITER_REF = 40000         # iteration cap of the reference runs (res_ctrl 1e-7)
# stalled: no new minimum of the largest residual for this many iterations (the environment
# variable exists for the local tests with the fake OpenFOAM only)
STALL_ITER = int(os.environ.get("AIRFOIL_STALL_ITER", 3000))
# A5: a settled run stops when its largest residual is at its target and the relative spread
# of L/D over the preceding SETTLE_WINDOW iterations is at most SETTLE_TOL (the environment
# variable exists for the local tests with the fake OpenFOAM only)
SETTLE_TOL = 1e-4
SETTLE_WINDOW = int(os.environ.get("AIRFOIL_SETTLE_WINDOW", 2000))
# A4: a run with n_surf above the Stage-1 list starts from the same design's solution on
# the SEQUENCE mesh of its own radius (mapFields); the environment variable exists for the
# local tests with the fake OpenFOAM only
SEQUENCE_ABOVE = int(os.environ.get("AIRFOIL_SEQUENCE_ABOVE", 256))
SEQUENCE = dict(n_surf=91, dy1=2e-5, res_ctrl=1e-6, gamg_tol=1e-8)
WRITE_EVERY = 1000           # checkpoint interval (iterations); purgeWrite keeps the latest
DIAG_WRITE_EVERY = 100       # numerics sweep: write interval; purgeWrite keeps the two latest
TARGETS = [0.01, 0.003]      # relative targets in L/D: primary 1 %, secondary 0.3 %
FALLBACK_TARGET = 0.03       # graded only if the 1 % target is not admitted (PROTOCOL.md 6)


@dataclass(frozen=True)
class Run:
    """One simulation: a design, a setting, the cores it runs on."""
    m: float
    t: float
    n_surf: int
    dy1: float
    radius: float
    res_ctrl: float
    gamg_tol: float
    np: int = 1
    numerics: str = ""          # "" : the configured variant (config/numerics.json)
    settle: bool = False        # A5: stops only when the residual target is met and L/D settled

    @property
    def num(self) -> str:
        """The numerics variant of this run (PROTOCOL.md 3.3)."""
        return self.numerics or numerics()

    @property
    def tag(self) -> str:
        """Deterministic, file-safe name.  Runs on more than one core carry _np<N> (their
        answer differs from the serial one only by round-off, but their timing is what
        Stage 0 measures); runs with the fallback numerics (PROTOCOL.md 3.3) carry _fb, runs
        with a variant of the numerics sweep _nv<variant>."""
        base = (f"m{self.m * 100:.1f}_t{self.t * 100:.0f}_n{self.n_surf}_y{self.dy1:.0e}"
                f"_R{self.radius:.0f}_rc{self.res_ctrl:.0e}_gt{self.gamg_tol:.0e}").replace("+", "")
        suffix = {"primary": "", "fallback": "_fb"}.get(self.num, f"_nv{self.num}")
        return base + (f"_np{self.np}" if self.np > 1 else "") + suffix + ("_st" if self.settle else "")

    @property
    def max_iter(self) -> int:
        return MAX_ITER_REF if self.res_ctrl < min(LADDERS[3]) else MAX_ITER

    def setting(self) -> dict:
        return {n: getattr(self, n) for n in NAMES}

    def to_dict(self) -> dict:
        d = asdict(self)
        d["tag"] = self.tag
        return d


def sequence_run(run):
    """The first run of a run finer than the Stage-1 lists (A4), or None; settled (A5)."""
    if run.n_surf <= SEQUENCE_ABOVE:
        return None
    return run_at((run.m, run.t), dict(SEQUENCE, radius=run.radius), numerics=run.numerics, settle=True)


def reference_run(design, setting=None):
    """A run of the reference setting (d), or of one of its siblings: settled (A5)."""
    return run_at(design, setting or reference_setting(), settle=True)


def far_field_run(design):
    """A5, item 3: the second run of the far-field correction, the first run of the design's
    reference (n_surf 91, dy1 2e-5) at R = 500 instead of 5000; settled."""
    return run_at(design, dict(SEQUENCE, radius=FAR_FIELD_RADIUS), settle=True)


FAR_FIELD_RADIUS = 500.0
FAR_FIELD_UNCERTAINTY = 0.2      # A5, item 4: share of the correction counted in u_adm


def far_field_correction(ld_5000, ld_500):
    """A5, item 3: Delta_R = (L/D(5000) - L/D(500)) / 9, first order in 1/R to R -> infinity."""
    return (ld_5000 - ld_500) * (1.0 / REFERENCE["radius"]) / (1.0 / FAR_FIELD_RADIUS - 1.0 / REFERENCE["radius"])


def run_at(design, setting, np=1, numerics="", settle=False) -> Run:
    m, t = design
    s = dict(setting)
    return Run(m=float(m), t=float(t), n_surf=int(s["n_surf"]), dy1=float(s["dy1"]),
               radius=float(s["radius"]), res_ctrl=float(s["res_ctrl"]),
               gamg_tol=float(s["gamg_tol"]), np=np, numerics=numerics, settle=settle)


def setting_of(idx) -> dict:
    return {n: LADDERS[k][j] for k, (n, j) in enumerate(zip(NAMES, idx))}


# ---------------------------------------------------------------- configuration files
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.environ.get("AIRFOIL_CONFIG", os.path.join(HERE, "config"))


def load_config(name):
    path = os.path.join(CONFIG_DIR, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def numerics():
    """'primary', or 'fallback' once the Stage-0 analysis has triggered it (config/numerics.json)."""
    c = load_config("numerics.json")
    return (c or {}).get("numerics") or "primary"


def hardest_design():
    """The design of Stage 1, from config/stage1.json (written by analysis/stage0.py)."""
    c = load_config("stage1.json")
    if not c or c.get("hardest") is None:
        raise RuntimeError("config/stage1.json has no hardest design yet: run the Stage-0 "
                           "analysis first (RUNBOOK.md, step 6)")
    return tuple(c["hardest"])


def radius_ladder():
    """The far-field list and the reference radius (the check of PROTOCOL.md 3.2 that could
    multiply them by ten is withdrawn by A5)."""
    return list(LADDERS[2]), REFERENCE["radius"]


def ladders_now():
    lad = [list(v) for v in LADDERS]
    lad[2] = radius_ladder()[0]
    return lad


def reference_setting():
    """(d): one step beyond the finest value of every list (PROTOCOL.md 5.1)."""
    return dict(REFERENCE, radius=radius_ladder()[1])


def default_start():
    """The user's default start: the middle of each list."""
    return {n: v[i] for n, v, i in zip(NAMES, ladders_now(), START_IDX)}


def tutorial_setting():
    """(b): the airFoil2D tolerances on the middle mesh (middle of the discretisation lists)."""
    return dict(default_start(), **TUTORIAL_TOLERANCES)


def cautious_setting():
    """(c): the finest value of every list."""
    return {n: v[-1] for n, v in zip(NAMES, ladders_now())}


# ---------------------------------------------------------------- sweeps
def sweep_check():
    """Two tiny runs for the 15-minute check job on the test partition."""
    s = dict(MID, n_surf=64, radius=50.0, res_ctrl=1e-3)
    return [run_at(CENTRE, s), run_at(CORNERS[-1], s)]


def sweep_meshcheck():
    """Meshes checkMesh must pass in the check job, before any other run (PROTOCOL.md 14, A2):
    the coarsest and the finest of the Stage-1 lists for the two most different corners, the
    reference, and the mesh-order check at n_surf 512."""
    lad = ladders_now()
    coarsest = {n: v[0] for n, v in zip(NAMES, lad)}
    finest = {n: v[-1] for n, v in zip(NAMES, lad)}
    ref = reference_setting()
    thick, thin = CORNERS[-1], CORNERS[0]
    return [run_at(thick, coarsest), run_at(thin, coarsest),
            run_at(thick, finest), run_at(thin, finest),
            run_at(thick, ref), run_at(thick, dict(ref, n_surf=L6["n_surf"]))]


def sweep_diag():
    """Diagnostic before Stage 0 (10 Nov 2026): the m 6 %, t 16 % corner oscillated with the
    primary numerics on the coarsest mesh.  Is it the numerics or the flow?  The corner on
    the coarsest mesh with the fallback numerics; the corner, and the thin cambered corner,
    at the middle setting with each variant (the primary runs at the middle setting are
    Stage-0 runs C and are reused there)."""
    coarse = dict(MID, n_surf=64, radius=50.0, res_ctrl=1e-3)
    return [run_at((0.06, 0.16), coarse, numerics="fallback"),
            run_at((0.06, 0.16), MID, numerics="primary"),
            run_at((0.06, 0.16), MID, numerics="fallback"),
            run_at((0.06, 0.08), MID, numerics="primary"),
            run_at((0.06, 0.08), MID, numerics="fallback")]


def sweep_stage0():
    """Stage 0 (PROTOCOL.md 4.1): cost, convergence, hardest design, far field."""
    runs = []
    # A: the centre design at three mesh levels, middle setting otherwise
    for n in (64, 128, 256):
        runs.append(run_at(CENTRE, dict(MID, n_surf=n)))
    # A': the centre design on the finest mesh at the finest tolerances
    runs.append(run_at(CENTRE, dict(MID, n_surf=256, res_ctrl=1e-6, gamg_tol=1e-8)))
    # B: cores -- the centre design at level 3, serial and on four cores
    runs.append(run_at(CENTRE, dict(MID, n_surf=181)))
    runs.append(run_at(CENTRE, dict(MID, n_surf=181), np=4))
    # C: the four corners at the middle setting and one level finer
    for d in CORNERS:
        for n in (128, 181):
            runs.append(run_at(d, dict(MID, n_surf=n)))
    # D: far field at the corners, level 1, finest tolerances: R = 50, 500, 5000 (A5)
    for d in CORNERS:
        for R in (50.0, 500.0, 5000.0):                    # A5: no 50,000; settled
            runs.append(run_at(d, dict(MID, n_surf=91, radius=R, res_ctrl=1e-6, gamg_tol=1e-8), settle=True))
    return runs


def sweep_stage1():
    """Stage 1: the 540-setting grid at the hardest design, ordered by mesh level
    (slowest index) so that each level is one contiguous range of array indices."""
    d = hardest_design()
    lad = ladders_now()
    return [run_at(d, {n: lad[k][j] for k, (n, j) in enumerate(zip(NAMES, idx))}) for idx in GRID]


def sweep_stage1ref():
    """The reference at the hardest design, its one-step-back siblings (PROTOCOL.md 5.2)
    and the mesh-order check at n_surf 512."""
    return stage1ref_runs(hardest_design())


def stage1ref_runs(d):
    """The runs of PROTOCOL.md 5.2 at design d (A5: settled, with the far-field run)."""
    ref = reference_setting()
    lad = ladders_now()
    runs = [reference_run(d, ref)]
    for name in ("n_surf", "dy1", "radius"):
        k = NAMES.index(name)
        # one step back: the finest list value, or the next one down where the reference
        # already sits at the finest (the radius, A5)
        back = lad[k][-1] if lad[k][-1] != ref[name] else lad[k][-2]
        runs.append(reference_run(d, dict(ref, **{name: back})))
    runs.append(reference_run(d, dict(ref, n_surf=L6["n_surf"])))
    runs.append(far_field_run(d))                       # A5: the far-field correction
    return runs


def stage2_reference_runs(s=None):
    """Stage 2, setting (d): the reference of every design and its far-field run (A5)."""
    return [reference_run(d, s) for d in DESIGNS] + [far_field_run(d) for d in DESIGNS]


def stage2_settings():
    c = load_config("stage2_settings.json")
    if not c or not c.get("settings"):
        raise RuntimeError("config/stage2_settings.json is not written yet: run the Stage-1 "
                           "analysis first (RUNBOOK.md, step 9)")
    return c["settings"]


def sweep_stage2():
    """Stage 2: the 25 designs at each setting of config/stage2_settings.json, setting by
    setting (a, b, c, d) and design by design."""
    runs = []
    for key in sorted(stage2_settings()):
        s = stage2_settings()[key]
        if key == "d":
            runs.extend(stage2_reference_runs(s))

        else:
            runs.extend(run_at(d, s) for d in DESIGNS)
    return runs


# the numerics sweep (10 Nov 2026, before Stage 0): on both meshes the primary variant
# never converges (largest residual 1e-3 to 1e-2, lift oscillating), the fallback diverges
# within 50 iterations, and the upper surface separates at 75-80 % chord at 4 degrees,
# with lift 25-40 % below the inviscid value.  Which numerics converge, and is the
# separation the mesh, the numerics or the turbulence model?
NUMERICS_DESIGNS = [(0.03, 0.12), (0.06, 0.16)]
VALIDATION = (0.0, 0.12)              # NACA 0012 with the closed trailing edge
NUMERICS_SETTING = dict(n_surf=64, dy1=2e-5, radius=50.0, res_ctrl=1e-6, gamg_tol=1e-8)


def sweep_numerics():
    """(a) every variant, one change at a time from the primary, on the coarsest mesh, at
    the centre design and the m 6 %, t 16 % corner; (b) mesh resolution n_surf 64, 128,
    256 with the primary and the robust variant, and NACA 0012 as a check against the
    known lift; (c) the Spalart-Allmaras model with the robust schemes, as a check of the
    k-omega SST set-up.  Far field R = 50, residualControl 1e-6, GAMG 1e-8, so that a
    variant that converges stops there and one that does not stalls."""
    s64 = dict(NUMERICS_SETTING)
    runs = []
    for v in ("primary", "limU", "upKW", "simple", "limU_upKW", "robust", "robust_nonorth", "upwind1"):
        for d in NUMERICS_DESIGNS:
            runs.append(run_at(d, s64, numerics=v))
    for n in (128, 256):
        for v in ("primary", "robust"):
            for d in NUMERICS_DESIGNS:
                runs.append(run_at(d, dict(s64, n_surf=n), numerics=v))
    for n in (64, 128, 256):
        runs.append(run_at(VALIDATION, dict(s64, n_surf=n), numerics="robust"))
    for n in (64, 128):
        for d in NUMERICS_DESIGNS + [VALIDATION]:
            runs.append(run_at(d, dict(s64, n_surf=n), numerics="sa"))
    return runs


# the convergence map (PROTOCOL.md 14, A3, "Before Stage 0"): with the A3 free stream and
# numerics the thick cambered sections converge on some meshes and level off at 2e-6 to
# 7e-5 on others.  Which of the 45 Stage-1 meshes converge at the tightest tolerances, and
# does the reference reach 1e-7?
CONVMAP_DESIGNS = [(0.06, 0.16), (0.06, 0.14), (0.045, 0.16), (0.045, 0.14)]


def sweep_convmap():
    """The 45 Stage-1 meshes (n_surf x dy1 x radius) of each design of CONVMAP_DESIGNS at
    residualControl 1e-6 and GAMG 1e-8, design by design; then the reference setting at
    m 6 %, t 16 % and m 4.5 %, t 16 %."""
    tight = dict(res_ctrl=min(LADDERS[3]), gamg_tol=min(LADDERS[4]))
    radii = radius_ladder()[0]
    runs = []
    for d in CONVMAP_DESIGNS:
        for n in LADDERS[0]:
            for dy in LADDERS[1]:
                for R in radii:
                    runs.append(run_at(d, dict(tight, n_surf=n, dy1=dy, radius=R)))
    for d in ((0.06, 0.16), (0.045, 0.16)):
        runs.append(run_at(d, reference_setting()))
    return runs


# the confirmation map (PROTOCOL.md 14, A4): the two cambered corners of the A4 design
# space on the 45 Stage-1 meshes at the tightest tolerances, and the reference at the
# thicker one, with the far-wake mesh of A4
CONFIRM_DESIGNS = [(0.04, 0.12), (0.04, 0.08)]


def sweep_confirm():
    """The 45 Stage-1 meshes of each design of CONFIRM_DESIGNS at residualControl 1e-6 and
    GAMG 1e-8, design by design; then the reference setting at m 4 %, t 12 %."""
    tight = dict(res_ctrl=min(LADDERS[3]), gamg_tol=min(LADDERS[4]))
    radii = radius_ladder()[0]
    runs = [run_at(d, dict(tight, n_surf=n, dy1=dy, radius=R))
            for d in CONFIRM_DESIGNS for n in LADDERS[0] for dy in LADDERS[1] for R in radii]
    runs.append(run_at(CONFIRM_DESIGNS[0], reference_setting()))
    return runs


# the settling measurement (12 November, after the confirmation map; no registered rule
# changes): the reference of m 4 %, t 12 % stalled with L/D in a slow damped oscillation of
# about +-2 %.  How many iterations does L/D need to settle, at R = 50,000 and at R = 5,000?
# Two runs of the reference mesh, started as in A4, that ignore the residual stop (a
# residualControl no run meets) and the stall watchdog, to the reference cap of 40,000
# iterations, with their whole force history kept.
SETTLE_DESIGN = (0.04, 0.12)
SETTLE_RES_CTRL = 1e-12
SETTLE_RADII = (50000.0, 5000.0)
# for the local tests with the fake OpenFOAM only: "n_surf,dy1" of a smaller mesh
_SETTLE_TEST = os.environ.get("AIRFOIL_SETTLE_TEST_MESH")


def sweep_settle():
    """The reference mesh at R = 50,000 and 5,000, to 40,000 iterations (no residual stop)."""
    mesh = dict(n_surf=REFERENCE["n_surf"], dy1=REFERENCE["dy1"])
    if _SETTLE_TEST:
        n, dy = _SETTLE_TEST.split(",")
        mesh = dict(n_surf=int(n), dy1=float(dy))
    return [run_at(SETTLE_DESIGN, dict(REFERENCE, radius=R, res_ctrl=SETTLE_RES_CTRL, **mesh))
            for R in SETTLE_RADII]


SWEEPS = {"check": sweep_check, "settle": sweep_settle, "numerics": sweep_numerics, "convmap": sweep_convmap, "confirm": sweep_confirm, "meshcheck": sweep_meshcheck, "diag": sweep_diag, "stage0": sweep_stage0, "stage1": sweep_stage1,
          "stage1ref": sweep_stage1ref, "stage2": sweep_stage2}


def sweep(name):
    if name not in SWEEPS:
        raise KeyError(f"unknown sweep {name!r}; known: {', '.join(SWEEPS)}")
    return SWEEPS[name]()


# ---------------------------------------------------------------- cost model
# Unmeasured until Stage 0 (COST.md): core-seconds per cell and iteration, and
# iterations to meet residualControl on the middle mesh, growing per level.
T_CELL_ITER = 4.6e-6         # measured: median of the convergence map (A4); was 3.0e-6
# A4: the growth per level as measured in the convergence map (1.29); A5: 1e-7 calibrated on
# the settled reference of the settling measurement (about 11,000 iterations at n_surf 362,
# three levels above 128, at R = 5000); the rest as first assumed.  Stage 0 recalibrates it.
ITERS_MID = {1e-3: 300, 1e-4: 700, 1e-5: 1200, 1e-6: 1800, 1e-7: 5100}
ITER_GROWTH = 1.29           # per mesh level (factor sqrt(2) in n_surf)


def predicted_iters(run: Run, iters_mid=None, growth=None) -> float:
    im = iters_mid or ITERS_MID
    g = growth or ITER_GROWTH
    level = math.log(run.n_surf / 128.0) / math.log(math.sqrt(2.0))
    base = im.get(run.res_ctrl)
    if base is None:
        keys = sorted(im)
        base = float(np.interp(math.log10(run.res_ctrl), [math.log10(k) for k in keys],
                               [im[k] for k in keys]))
    return min(base * g ** level, run.max_iter)
