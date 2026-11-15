#!/usr/bin/env python3
"""A stand-in for the OpenFOAM executables the driver calls, for local tests only.

Invoked through symlinks named checkMesh, simpleFoam, decomposePar,
reconstructPar and mpirun.  It writes what the driver reads, in OpenFOAM's
formats: a solver log with residuals and ExecutionTime, forceCoeffs output
(coefficient.dat), time directories honouring writeInterval / purgeWrite /
startFrom latestTime / stopAt writeNow, and, in -postProcess mode, the yPlus,
wallShearStress and C fields with the aerofoil patch values.

Behaviour switches (environment): FAKEFOAM_STALL=1 (residuals plateau),
FAKEFOAM_CRASH_AT=k (floating-point exception at iteration k),
FAKEFOAM_SETUP_FAIL=1 (fatal error before the first iteration),
FAKEFOAM_SLEEP=s (seconds per iteration).
"""
import math
import os
import re
import shutil
import sys
import time

prog = os.path.basename(sys.argv[0])
case = os.getcwd()
BANNER = """/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Version:  v2606
\\*---------------------------------------------------------------------------*/
Build  : fake
"""


def rd(rel):
    return open(os.path.join(case, rel)).read()


def entry(txt, key):
    m = re.search(r"^\s*%s\s+([^;]+);" % re.escape(key), txt, flags=re.M)
    return m.group(1).strip() if m else None


def cells_from_mesh():
    m = re.search(r"nCells:(\d+)", open(os.path.join(case, "constant/polyMesh/owner")).read(400))
    return int(m.group(1))


def n_surf():
    b = rd("constant/polyMesh/boundary")
    return int(re.search(r"airfoil\s*\{[^}]*nFaces\s+(\d+);", b).group(1)) // 2


def times(d):
    return sorted((float(x), x) for x in os.listdir(d) if re.fullmatch(r"[0-9.eE+-]+", x)
                  and os.path.isdir(os.path.join(d, x)))


def check_mesh():
    print(BANNER)
    n = cells_from_mesh()
    print(f"Mesh stats\n    points:           {2 * n}\n    faces:            {4 * n}\n    cells:            {n}\n")
    # the lines of a v2606 checkMesh, with their full stops after numbers
    print("Checking geometry...\n    Overall domain bounding box (-49.5 -50 0) (50.5 50 1)")
    print(" ***High aspect ratio cells found, Max aspect ratio: 81062.3, number of cells 12")
    print("  <<Writing 12 cells with high aspect ratio to set highAspectRatioCells")
    print("    Minimum face area = 1.2e-09. Maximum face area = 54.1.  Face area magnitudes OK.")
    print("    Min volume = 4.17398418547e-08. Max volume = 54.1.  Total volume = 9999.9.  Cell volumes OK.")
    print("    Mesh non-orthogonality Max: 21.3 average: 4.1\n    Non-orthogonality check OK.")
    print("    Face pyramids OK.\n    Max skewness = 0.61 OK.\n    Coupled point location match (average 0) OK.")
    if os.environ.get("FAKEFOAM_CHECKMESH_FAIL"):
        print(" ***Error in face pyramids: 4 faces are incorrectly oriented.")
        print("\nFailed 2 mesh checks.\n\nEnd")
        return
    print("\nFailed 1 mesh checks.\n\nEnd")


def write_fields(tdir, it, cells):
    """Ascii fields as simpleFoam writes them: a nonuniform internalField per field; the
    values drift with the iteration so that two written times differ."""
    os.makedirs(tdir, exist_ok=True)
    names = ("U", "p", "k", "omega", "nut")
    if os.path.exists(os.path.join(case, "0", "nuTilda")):
        names += ("nuTilda",)
    for f in names:
        vec = f == "U"
        if vec:
            body = "".join(f"({1.0 + 1e-3 * math.sin(0.01 * it + 1e-3 * i):.9g} {0.02 * math.cos(1e-3 * i):.6g} 0)\n"
                           for i in range(cells))
        else:
            body = "".join(f"{1e-4 * (1 + math.sin(1e-3 * i)) + 1e-6 * it:.9g}\n" for i in range(cells))
        kind = "vector" if vec else "scalar"
        with open(os.path.join(tdir, f), "w") as fh:
            fh.write(f"FoamFile\n{{\n    class       vol{kind.capitalize()}Field;\n    object      {f};\n}}\n"
                     f"dimensions      [0 0 0 0 0 0 0];\n\n"
                     f"internalField   nonuniform List<{kind}> \n{cells}\n(\n{body})\n;\n\n"
                     "boundaryField\n{\n    airfoil { type calculated; value uniform 0; }\n}\n")


def simple_foam():
    print(BANNER)
    if os.environ.get("FAKEFOAM_SETUP_FAIL"):
        print("--> FOAM FATAL IO ERROR: keyword decayControl is undefined")
        sys.exit(1)
    print("Selecting RAS turbulence model kOmegaSST\nkOmegaSSTCoeffs\n{\n    alphaK1         0.85;\n"
          "    decayControl    true;\n    kInf            4e-07;\n    omegaInf        133.333;\n}\n")
    parallel = "-parallel" in sys.argv
    root = os.path.join(case, "processor0") if parallel else case
    if parallel and not os.path.isdir(root):
        os.makedirs(os.path.join(root, "0"))
    ctl = rd("system/controlDict")
    end = int(entry(ctl, "endTime"))
    every = int(entry(ctl, "writeInterval"))
    purge = int(entry(ctl, "purgeWrite") or 0)
    sol = rd("system/fvSolution")
    rc = float(re.search(r"residualControl\s*\{\s*p\s+([0-9.eE+-]+);", sol).group(1))
    ts = times(root)
    it = int(ts[-1][0]) if ts else 0
    ns = n_surf()
    cells = cells_from_mesh()
    # the answer: converges with the mesh (second order) and with the residual level
    cl_inf, cd_inf = 0.95, 0.0090
    cl_t = cl_inf * (1 - 0.8 * (64.0 / ns) ** 2 * 0.05)
    cd_t = cd_inf * (1 + 0.6 * (64.0 / ns) ** 2 * 0.05)
    pp = os.path.join(case, "postProcessing/forceCoeffs", str(it))
    os.makedirs(pp, exist_ok=True)
    fc = open(os.path.join(pp, "coefficient.dat"), "a")
    fc.write("# Force and moment coefficients\n# dragDir : (1 0 0)\n#\n"
             "# Time          \tCd\tCd(f)\tCd(r)\tCl\tCl(f)\tCl(r)\tCmPitch\tCmRoll\tCmYaw\tCs\tCs(f)\tCs(r)\n")
    crash = int(os.environ.get("FAKEFOAM_CRASH_AT", "0"))
    stall = bool(os.environ.get("FAKEFOAM_STALL"))
    sleep = float(os.environ.get("FAKEFOAM_SLEEP", "0"))
    exe = 0.0
    while it < end:
        it += 1
        r = 0.2 * 0.985 ** it
        if stall:
            r = max(r, 3e-4 * (1 + 0.1 * math.sin(it)))
        print(f"Time = {it}\n")
        for fld, f in (("Ux", 1.0), ("Uy", 2.0), ("p", 3.0), ("omega", 0.5), ("k", 0.7)):
            solver = "GAMG" if fld == "p" else "smoothSolver"
            print(f"{solver}:  Solving for {fld}, Initial residual = {r * f / 3:.6g}, "
                  f"Final residual = {r * f / 30:.6g}, No Iterations 2")
        exe += cells * 3e-6
        print(f"ExecutionTime = {exe:.2f} s  ClockTime = {int(exe) + 1} s\n")
        cd = cd_t * (1 + 0.3 * r)
        cl = cl_t * (1 - 0.3 * r)
        fc.write(f"{it}\t{cd:.10g}\t0\t0\t{cl:.10g}\t0\t0\t-0.1\t0\t0\t0\t0\t0\n")
        fc.flush()
        sys.stdout.flush()
        if crash and it >= crash:
            print("#0  Foam::error::printStack(Foam::Ostream&)\n#1  Foam::sigFpe::sigHandler(int)")
            print("Floating point exception (core dumped)")
            sys.exit(136)
        if sleep:
            time.sleep(sleep)
        converged = r * 3 / 3 < rc and r * 2 / 3 < rc and r / 3 < rc
        stop_now = entry(rd("system/controlDict"), "stopAt") == "writeNow"
        if stop_now and os.environ.get("FAKEFOAM_IO_FAIL_ON_STOP"):
            # what v2606 did on Hamilton when controlDict was caught half-written
            print("--> FOAM FATAL IO ERROR: (openfoam-2606)\nproblem while reading header for object controlDict")
            print("\nFOAM exiting\n")
            sys.exit(1)
        if converged or stop_now or it % every == 0 or it == end:
            write_fields(os.path.join(root, str(it)), it, cells)
            if purge:
                for t, name in [x for x in times(root) if x[0] > 0][:-purge]:
                    shutil.rmtree(os.path.join(root, name))     # purgeWrite: keep the latest

        if converged:
            print(f"\nSIMPLE solution converged in {it} iterations\n")
            break
        if stop_now:
            break
    print("End")


def post_process():
    print(BANNER)
    t = times(case)[-1][1]
    ns = n_surf()
    nf = 2 * ns
    th = [2 * math.pi * (i + 0.5) / nf for i in range(nf)]
    x = [0.5 + 0.5 * math.cos(a) for a in th]
    y = [0.06 * math.sin(a) for a in th]
    yp = [0.4 + 0.5 * abs(math.sin(a)) for a in th]
    tw = [(0.0016 * (1 + abs(math.cos(a))), 0.0, 0.0) for a in th]

    def wr(name, cls, body):
        with open(os.path.join(case, t, name), "w") as fh:
            fh.write("FoamFile { version 2.0; format ascii; class %s; object %s; }\n" % (cls, name))
            fh.write("dimensions [0 0 0 0 0 0 0];\ninternalField uniform 0;\nboundaryField\n{\n")
            fh.write(body + "    frontAndBack { type empty; }\n}\n")

    def lst(vals, vec=False):
        if vec:
            return "nonuniform List<vector> %d\n(\n%s\n)\n;" % (len(vals), "\n".join("(%r %r %r)" % v for v in vals))
        return "nonuniform List<scalar> %d\n(\n%s\n)\n;" % (len(vals), "\n".join(repr(v) for v in vals))

    wr("yPlus", "volScalarField", "    farfield { type calculated; value uniform 0; }\n"
       "    airfoil { type calculated; value %s }\n" % lst(yp))
    wr("wallShearStress", "volVectorField", "    farfield { type calculated; value uniform (0 0 0); }\n"
       "    airfoil { type calculated; value %s }\n" % lst(tw, True))
    wr("C", "volVectorField", "    airfoil { type calculated; value %s }\n"
       "    farfield { type calculated; value uniform (0 0 0); }\n" % lst([(a, b, 0.5) for a, b in zip(x, y)], True))
    print("End")


def reconstruct():
    src = os.path.join(case, "processor0")
    t = times(src)[-1][1]
    shutil.copytree(os.path.join(src, t), os.path.join(case, t), dirs_exist_ok=True)
    print("End")


if prog == "checkMesh":
    check_mesh()
elif prog == "simpleFoam":
    post_process() if "-postProcess" in sys.argv else simple_foam()
elif prog == "mapFields":
    # A4: map the solution of a coarser case (argument 1) onto this one; needs the dictionary
    print(BANNER)
    if not os.path.exists(os.path.join(case, "system", "mapFieldsDict")):
        print("--> FOAM FATAL IO ERROR: cannot find file system/mapFieldsDict")
        sys.exit(1)
    src = os.path.join(case, sys.argv[1])
    if not times(src) or times(src)[-1][0] <= 0:
        print(f"--> FOAM FATAL ERROR: no time to map in {src}")
        sys.exit(1)
    print(f"Mapping fields for time {times(src)[-1][1]}\n\n    interpolating U\n\nEnd")
elif prog == "decomposePar":
    os.makedirs(os.path.join(case, "processor0", "0"), exist_ok=True)
    print("End")
elif prog == "reconstructPar":
    reconstruct()
elif prog == "mpirun":
    args = [a for a in sys.argv[1:]]
    i = args.index("simpleFoam")
    os.execv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "simpleFoam"), args[i:])
else:
    sys.exit(f"fakefoam: unknown program {prog}")
