"""The OpenFOAM case of one run: simpleFoam, k-omega SST, wall-resolved.

Steady incompressible RANS at Re = 3e6 on chord 1, free stream (1, 0, 0) with the
section pitched by alpha.  Boundary conditions: noSlip on the aerofoil with the
low-Re wall treatment (nut = 0 through nutLowReWallFunction, k = 0, omega from
omegaWallFunction, which takes its viscous-sublayer value below y+ = 11);
freestreamVelocity / freestreamPressure on the whole far field; inletOutlet for k
and omega.  Free-stream turbulence (TMR's SST-sust values, PROTOCOL.md 14, A3) is
sustained by the decay control of kOmegaSST, so it reaches the aerofoil unchanged
whatever the far-field radius.

The numerics (PROTOCOL.md 2, as amended by A3): "primary", SIMPLEC with second-order
upwind for U and first-order upwind for k and omega.  The registered "fallback" (SIMPLE
with the airFoil2D relaxation) is withdrawn by A3; the other variants serve the numerics
sweep of 10 November.
"""
from __future__ import annotations

import os

from . import study
from .cgrid import SPAN, CGrid

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  errbudget airfoil study (PROTOCOL.md); written by af/foamcase.py
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def _write(case, rel, cls, obj, body):
    """Write a dictionary atomically: OpenFOAM re-reads controlDict while it runs
    (runTimeModifiable), and a file caught half-written stops it with a fatal I/O error."""
    path = os.path.join(case, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = os.path.join(os.path.dirname(path), f".{os.path.basename(path)}.tmp")
    with open(tmp, "w") as fh:
        fh.write(HEADER.format(cls=cls, obj=obj) + body)
    os.replace(tmp, path)


def control_dict(run, start="latestTime", stop="endTime", cofr=(0.25, 0.0), diag=False):
    """diag: the numerics sweep writes every DIAG_WRITE_EVERY iterations and keeps the two
    latest times, from which the driver measures where the solution still changes."""
    every, purge = (study.DIAG_WRITE_EVERY, 2) if diag else (study.WRITE_EVERY, 1)
    return f"""
application     simpleFoam;
startFrom       {start};
startTime       0;
stopAt          {stop};
endTime         {run.max_iter};
deltaT          1;
writeControl    timeStep;
writeInterval   {every};
purgeWrite      {purge};
writeFormat     ascii;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;

functions
{{
    forceCoeffs
    {{
        type            forceCoeffs;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   1;
        log             false;
        patches         (airfoil);
        rho             rhoInf;
        rhoInf          1;
        CofR            ({float(cofr[0])!r} {float(cofr[1])!r} {SPAN / 2!r});
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        magUInf         {study.U_INF!r};
        lRef            {study.CHORD!r};
        Aref            {study.CHORD * SPAN!r};
    }}
    // the three below run at write times and in `simpleFoam -postProcess -latestTime`,
    // which writes the fields yPlus, wallShearStress and C read by the driver
    yPlusRun
    {{
        type            yPlus;
        libs            (fieldFunctionObjects);
        executeControl  writeTime;
        writeControl    writeTime;
        log             false;
    }}
    wallShearStressRun
    {{
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        patches         (airfoil);
        executeControl  writeTime;
        writeControl    writeTime;
        log             false;
    }}
    cellCentresRun
    {{
        type            writeCellCentres;
        libs            (fieldFunctionObjects);
        executeControl  writeTime;
        writeControl    writeTime;
    }}
}}
"""


# The numerics variants (PROTOCOL.md 3.3).  "primary" and "fallback" are the registered
# pair; the others exist for the numerics sweep (study.sweep_numerics), which asks why the
# primary variant stalls and the fallback diverges.  Each entry: pressure-velocity
# coupling (SIMPLEC or SIMPLE), relaxation factors (p None: none), gradient of U and of
# k, omega (cellLimited or plain Gauss linear), convection of U and of k, omega
# (linearUpwind or upwind), non-orthogonal correctors.
VARIANTS = {
    # A3 (PROTOCOL.md 14): the primary variant with first-order convection of k and omega;
    # as registered it had linearUpwind ("primary_v1").  The fallback is the registered one,
    # withdrawn by A3 (it diverged in every run): kept for the record and the diag sweep.
    "primary":        dict(consistent=True, p=None, U=0.9, kw=0.7, grad_U="linear", grad_kw="limited",
                           div_U="linearUpwind", div_kw="upwind", nonorth=0),
    "fallback":       dict(consistent=False, p=0.3, U=0.7, kw=0.7, grad_U="linear", grad_kw="linear",
                           div_U="linearUpwind", div_kw="upwind", nonorth=0),
    "primary_v1":     dict(consistent=True, p=None, U=0.9, kw=0.7, grad_U="linear", grad_kw="limited",
                           div_U="linearUpwind", div_kw="linearUpwind", nonorth=0),
    # the numerics sweep of 10 November: one change at a time from primary_v1
    "limU":           dict(consistent=True, p=None, U=0.9, kw=0.7, grad_U="limited", grad_kw="limited",
                           div_U="linearUpwind", div_kw="linearUpwind", nonorth=0),
    "upKW":           dict(consistent=True, p=None, U=0.9, kw=0.7, grad_U="linear", grad_kw="limited",
                           div_U="linearUpwind", div_kw="upwind", nonorth=0),
    "simple":         dict(consistent=False, p=0.3, U=0.7, kw=0.7, grad_U="linear", grad_kw="limited",
                           div_U="linearUpwind", div_kw="linearUpwind", nonorth=0),
    # combinations
    "limU_upKW":      dict(consistent=True, p=None, U=0.9, kw=0.7, grad_U="limited", grad_kw="limited",
                           div_U="linearUpwind", div_kw="upwind", nonorth=0),
    "robust":         dict(consistent=True, p=None, U=0.7, kw=0.5, grad_U="limited", grad_kw="limited",
                           div_U="linearUpwind", div_kw="upwind", nonorth=0),
    "robust_nonorth": dict(consistent=True, p=None, U=0.7, kw=0.5, grad_U="limited", grad_kw="limited",
                           div_U="linearUpwind", div_kw="upwind", nonorth=1),
    # first order throughout: a diagnostic, never a candidate
    "upwind1":        dict(consistent=True, p=None, U=0.7, kw=0.5, grad_U="limited", grad_kw="limited",
                           div_U="upwind", div_kw="upwind", nonorth=0),
    # the Spalart-Allmaras model (screening of 11 November): the airFoil2D tutorial's
    # numerics (SIMPLE, p 0.3, U and nuTilda 0.7, linearUpwind U and nuTilda), and SIMPLEC
    # with the relaxation of the primary variant
    "sa_tut":         dict(consistent=False, p=0.3, U=0.7, kw=0.7, grad_U="linear", grad_kw="linear",
                           div_U="linearUpwind", div_kw="linearUpwind", nonorth=0, model="SpalartAllmaras"),
    "sa_c":           dict(consistent=True, p=None, U=0.9, kw=0.7, grad_U="linear", grad_kw="linear",
                           div_U="linearUpwind", div_kw="linearUpwind", nonorth=0, model="SpalartAllmaras"),
    # the Spalart-Allmaras model with the robust schemes: a check of the SST set-up
    "sa":             dict(consistent=True, p=None, U=0.7, kw=0.5, grad_U="limited", grad_kw="limited",
                           div_U="linearUpwind", div_kw="upwind", nonorth=0, model="SpalartAllmaras"),
}
NU_TILDA_INF = 3.0 * study.NU          # NASA TMR free-stream value for SA (nu_t / nu = 0.21)


def model_of(numerics):
    return VARIANTS[numerics].get("model", "kOmegaSST")


def fv_schemes(numerics="primary"):
    v = VARIANTS[numerics]
    if v["grad_U"] == "linear" and v["grad_kw"] == "linear":
        grad = "gradSchemes     { default Gauss linear; }\n"
    else:
        grad = "gradSchemes\n{\n    default         Gauss linear;\n"
        if v["grad_U"] == "limited":
            grad += "    grad(U)         cellLimited Gauss linear 1;\n"
        if v["grad_kw"] == "limited":
            grad += "    grad(k)         cellLimited Gauss linear 1;\n    grad(omega)     cellLimited Gauss linear 1;\n"
        grad += "}\n"

    def div(field, how):
        return f"bounded Gauss linearUpwind grad({field})" if how == "linearUpwind" else "bounded Gauss upwind"
    turb = ["k", "omega"] if model_of(numerics) == "kOmegaSST" else ["nuTilda"]
    if model_of(numerics) != "kOmegaSST" and v["grad_kw"] == "limited":
        grad = grad.replace("    grad(k)         cellLimited Gauss linear 1;\n    grad(omega)     cellLimited Gauss linear 1;\n",
                            "    grad(nuTilda)   cellLimited Gauss linear 1;\n")
    return ("\nddtSchemes      { default steadyState; }\n" + grad
            + "divSchemes\n{\n    default         none;\n"
            + f"    div(phi,U)      {div('U', v['div_U'])};\n"
            + "".join(f"    {'div(phi,' + f + ')':<15} {div(f, v['div_kw'])};\n" for f in turb)
            + "    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n"
            + "laplacianSchemes { default Gauss linear corrected; }\n"
            + "interpolationSchemes { default linear; }\n"
            + "snGradSchemes   { default corrected; }\n"
            + "wallDist        { method meshWave; }\n")


def fv_solution(run, numerics="primary"):
    rc, gt = run.res_ctrl, run.gamg_tol
    if getattr(run, "settle", False):
        rc = 1e-12          # A5: OpenFOAM's residual stop off; the driver stops a settled run
    v = VARIANTS[numerics]
    turb = "(k|omega)" if model_of(numerics) == "kOmegaSST" else "nuTilda"
    algo = f"    consistent      {'yes' if v['consistent'] else 'no'};\n"
    fields = "" if v["p"] is None else f"    fields\n    {{\n        p               {v['p']};\n    }}\n"
    relax = (f"relaxationFactors\n{{\n{fields}    equations\n    {{\n"
             f"        U               {v['U']};\n        {chr(34) + turb + chr(34):<16}{v['kw']};\n    }}\n}}\n")
    text = f"""
solvers
{{
    p
    {{
        solver          GAMG;
        smoother        GaussSeidel;
        tolerance       {gt!r};
        relTol          0.1;
    }}
    "(U|k|omega)"
    {{
        solver          smoothSolver;
        smoother        GaussSeidel;
        nSweeps         2;
        tolerance       {gt * study.SOLVER_RATIO!r};
        relTol          0.1;
    }}
}}

SIMPLE
{{
    nNonOrthogonalCorrectors {v['nonorth']};
{algo}    residualControl
    {{
        p               {rc!r};
        U               {rc!r};
        "(k|omega)"     {rc!r};
    }}
}}

{relax}"""
    if turb == "nuTilda":
        text = text.replace('"(U|k|omega)"', '"(U|nuTilda)"').replace('"(k|omega)"     ', '"nuTilda"       ')
    return text


def initial_fields():
    k, w = study.K_INF, study.OMEGA_INF
    U = f"({study.U_INF!r} 0 0)"
    return {
        "U": ("volVectorField", f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform {U};
boundaryField
{{
    airfoil         {{ type noSlip; }}
    farfield        {{ type freestreamVelocity; freestreamValue uniform {U}; value uniform {U}; }}
    frontAndBack    {{ type empty; }}
}}
"""),
        "p": ("volScalarField", """
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    airfoil         { type zeroGradient; }
    farfield        { type freestreamPressure; freestreamValue uniform 0; value uniform 0; }
    frontAndBack    { type empty; }
}
"""),
        "k": ("volScalarField", f"""
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform {k!r};
boundaryField
{{
    airfoil         {{ type fixedValue; value uniform 0; }}
    farfield        {{ type inletOutlet; inletValue uniform {k!r}; value uniform {k!r}; }}
    frontAndBack    {{ type empty; }}
}}
"""),
        "omega": ("volScalarField", f"""
dimensions      [0 0 -1 0 0 0 0];
internalField   uniform {w!r};
boundaryField
{{
    airfoil         {{ type omegaWallFunction; value uniform {w!r}; }}
    farfield        {{ type inletOutlet; inletValue uniform {w!r}; value uniform {w!r}; }}
    frontAndBack    {{ type empty; }}
}}
"""),
        "nuTilda": ("volScalarField", f"""
dimensions      [0 2 -1 0 0 0 0];
internalField   uniform {NU_TILDA_INF!r};
boundaryField
{{
    airfoil         {{ type fixedValue; value uniform 0; }}
    farfield        {{ type inletOutlet; inletValue uniform {NU_TILDA_INF!r}; value uniform {NU_TILDA_INF!r}; }}
    frontAndBack    {{ type empty; }}
}}
"""),
        "nut": ("volScalarField", """
dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;
boundaryField
{
    airfoil         { type nutLowReWallFunction; value uniform 0; }
    farfield        { type calculated; value uniform 0; }
    frontAndBack    { type empty; }
}
"""),
    }


def turbulence_properties(model="kOmegaSST"):
    if model != "kOmegaSST":
        return f"""
simulationType  RAS;

RAS
{{
    RASModel        {model};
    turbulence      on;
    printCoeffs     on;
}}
"""
    return f"""
simulationType  RAS;

RAS
{{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
    kOmegaSSTCoeffs
    {{
        // sustaining terms (Spalart & Rumsey 2007): free-stream k and omega do not decay
        decayControl    yes;
        kInf            {study.K_INF!r};
        omegaInf        {study.OMEGA_INF!r};
    }}
}}
"""


def write_case(case, run, grid: CGrid, numerics="primary", diag=False):
    """Write every input file of the case and the mesh (constant/polyMesh)."""
    qc = grid.foil.quarter_chord()
    grid.write_polymesh(case)
    _write(case, "system/controlDict", "dictionary", "controlDict", control_dict(run, cofr=qc, diag=diag))
    _write(case, "system/fvSchemes", "dictionary", "fvSchemes", fv_schemes(numerics))
    _write(case, "system/fvSolution", "dictionary", "fvSolution", fv_solution(run, numerics))
    _write(case, "system/decomposeParDict", "dictionary", "decomposeParDict",
           f"\nnumberOfSubdomains {max(run.np, 1)};\nmethod scotch;\n")
    _write(case, "constant/transportProperties", "dictionary", "transportProperties",
           f"\ntransportModel  Newtonian;\nnu              {study.NU!r};\n")
    _write(case, "constant/turbulenceProperties", "dictionary", "turbulenceProperties",
           turbulence_properties(model_of(numerics)))
    for name, (cls, body) in initial_fields().items():
        if name == "nuTilda" and model_of(numerics) == "kOmegaSST":
            continue                    # the SST cases are written as registered
        _write(case, f"0/{name}", cls, name, body)


def write_map_fields_dict(case):
    """system/mapFieldsDict for mapping a solution from another mesh (A4): internal fields
    only; the patches keep their own conditions."""
    _write(case, "system/mapFieldsDict", "dictionary", "mapFieldsDict",
           "\npatchMap        ( );\ncuttingPatches  ( );\n")


def set_control(case, run, stop="endTime", cofr=(0.25, 0.0), diag=False):
    """Rewrite controlDict (resume from the latest time; stopAt writeNow to stop cleanly)."""
    _write(case, "system/controlDict", "dictionary", "controlDict",
           control_dict(run, stop=stop, cofr=cofr, diag=diag))
