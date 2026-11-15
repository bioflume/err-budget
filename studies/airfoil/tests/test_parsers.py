"""Readers of OpenFOAM v2606 output, on text in its real format."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from af import foamio  # noqa: E402
from check_smoke import TUTORIAL, tutorial_values  # noqa: E402

CHECKMESH = """Mesh stats
    points:           27940
    faces:            55056
    internal faces:   27116
    cells:            13664
    boundary patches: 3
Checking geometry...
    Overall domain bounding box (-49.5 -50 0) (50.5 50 1)
    Mesh has 2 geometric (non-empty/wedge) directions (1 1 0)
    Max cell openness = 1.9e-16 OK.
 ***High aspect ratio cells found, Max aspect ratio: 1.04165e+06, number of cells 1302
  <<Writing 1302 cells with high aspect ratio to set highAspectRatioCells
    Minimum face area = 1.99e-05. Maximum face area = 54.1.  Face area magnitudes OK.
    Min volume = 4.17398418547e-08. Max volume = 54.1.  Total volume = 7853.9.  Cell volumes OK.
    Mesh non-orthogonality Max: 19.2479 average: 3.1
    Non-orthogonality check OK.
    Face pyramids OK.
    Max skewness = 0.583 OK.

Failed 1 mesh checks.
"""

# tutorials/incompressible/simpleFoam/airFoil2D/system/fvSolution of v2606, verbatim from 'solvers'
TUTORIAL_FILE = """FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.1;
        smoother        GaussSeidel;
    }

    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        nSweeps         2;
        tolerance       1e-08;
        relTol          0.1;
    }

    nuTilda
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        nSweeps         2;
        tolerance       1e-08;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;

    residualControl
    {
        p               1e-5;
        U               1e-5;
        nuTilda         1e-5;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        nuTilda         0.7;
    }
}
"""


def test_check_mesh_numbers_with_full_stops():
    c = foamio.parse_check_mesh(CHECKMESH)
    assert c["cells"] == 13664 and c["points"] == 27940
    assert c["min_volume"] == pytest.approx(4.17398418547e-08) and c["max_volume"] == pytest.approx(54.1)
    assert c["max_aspect"] == pytest.approx(1.04165e+06)
    assert c["max_non_orth"] == pytest.approx(19.2479) and c["avg_non_orth"] == pytest.approx(3.1)
    assert c["max_skewness"] == pytest.approx(0.583)
    assert c["failed_checks"] == 1 and not c["fatal"] and not c["mesh_ok"]
    assert foamio.parse_check_mesh("    Max aspect ratio = 812.5 OK.\n\nMesh OK.\n")["max_aspect"] == 812.5


@pytest.mark.parametrize("text", [TUTORIAL_FILE,
                                  TUTORIAL_FILE.replace("    U\n    {", '    "(U|nuTilda)"\n    {')])
def test_tutorial_file(tmp_path, text):
    p = tmp_path / "fvSolution"
    p.write_text(text)
    v = tutorial_values(str(p))
    assert v == TUTORIAL


@pytest.mark.parametrize("text,kind", [
    ("Starting time loop\nTime = 1\n#1  Foam::sigFpe::sigHandler(int)\n", "numerical"),
    ("Time = 3\n--> FOAM FATAL ERROR: Maximum number of iterations exceeded\n", "numerical"),
    ("Time = 6197\n--> FOAM FATAL IO ERROR: (openfoam-2606)\nproblem while reading header\n", "io"),
    ("Create mesh\n--> FOAM FATAL IO ERROR: keyword decayControl is undefined\n", "setup"),
    ("Time = 1\nExecutionTime = 0.1 s  ClockTime = 0 s\nEnd\n", None),
])
def test_failure_kinds(tmp_path, text, kind):
    p = tmp_path / "log"
    p.write_text(text)
    assert foamio.log_failure_kind(str(p)) == kind


def test_check_mesh_failures_allow_only_the_aspect_ratio():
    assert foamio.check_mesh_failures(CHECKMESH) == []
    bad = CHECKMESH.replace("    Face pyramids OK.", " ***Error in face pyramids: 4 faces are incorrectly oriented.")
    bad = bad.replace("Failed 1 mesh checks.", "Failed 2 mesh checks.")
    assert foamio.check_mesh_failures(bad) == [" ***Error in face pyramids: 4 faces are incorrectly oriented.".strip()]
    sev = CHECKMESH.replace("    Non-orthogonality check OK.", "   *Number of severely non-orthogonal (> 70 degrees) faces: 12.")
    assert len(foamio.check_mesh_failures(sev)) == 1
    assert foamio.check_mesh_failures("Failed 3 mesh checks.\n") == ["checkMesh: Failed 3 mesh checks"]
