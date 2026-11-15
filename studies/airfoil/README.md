# airfoil: the OpenFOAM use case of the errbudget article

errbudget tunes the numerical setting of OpenFOAM (`simpleFoam`, k-ω SST) once, at the
hardest of 25 NACA 4-digit sections, for L/D to 1 %; a 25-design camber × thickness
optimisation then reuses the setting. Graded on stored runs against a reference, and
compared with the tutorial's setting, a cautious one, and the procedures of practice.

| File | What |
|---|---|
| `PROTOCOL.md` | the pre-registered protocol (committed before any script) |
| `RUNBOOK.md` | the commands on Hamilton, in order, stage by stage |
| `COST.md` | the cost estimate and its assumptions |
| `af/naca.py` | NACA 4-digit geometry, closed trailing edge, pitched by α |
| `af/grading.py`, `af/cgrid.py` | point distributions; the C-grid computed point by point, its check as checkMesh measures, and the polyMesh writer |
| `af/foamcase.py`, `af/foamio.py` | the OpenFOAM case files; readers of logs, forces, patch fields |
| `af/study.py`, `af/costmodel.py` | parameters, settings, sweeps, tags; cells and predicted cost |
| `driver.py` | one run by sweep and index (`--sweep --index --list --dry-run`), or by setting (`--run-setting`) |
| `errbudget_job.yaml`, `live_tune.py` | errbudget's job file for this case, and a runner through its Python API |
| `slurm/` | configuration, environment, setup, check job, array job, submit, status, pack |
| `analysis/` | check, accounting, Stages 0–2, the replay copied from the article, the cost estimate, synthetic runs |
| `config/` | written by the analyses: the hardest design (Stage 0), the Stage-2 settings (Stage 1) |
| `tests/` | geometry; the grid and the polyMesh writer; the driver against a fake OpenFOAM; the analyses on synthetic runs |
| `notes/` | dated notes of diagnoses (e.g. the free stream and the numerics, 10 November) |
| `dist/` | `airfoil-batch.tar.gz` for Hamilton (`bash make_dist.sh`) and its SHA-256 |

Local checks (about a minute; NumPy, SciPy, Matplotlib, PyYAML and pytest; the analysis tests
use the engine in the repository's `src/`, or `$ERRBUDGET_SRC`):

```bash
cd airfoil/tests && PYTHONDONTWRITEBYTECODE=1 /opt/anaconda3/bin/python3 -m pytest -q -p no:cacheprovider
for f in ../slurm/*.sh ../slurm/*.slurm ../make_dist.sh; do bash -n "$f"; done
```

The article's repository documents an `errbudget` command line (`errbudget run job.yaml`),
but its `src/errbudget` has the engine and the Python API only. `live_tune.py` therefore
reads the job file and calls `errbudget.tune`; with a release that has the command line,
`errbudget run errbudget_job.yaml` should give the same result.
