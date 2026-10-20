# errbudget-paper: the article's code, data and reproduction

This repository reproduces the article

> G. Kabacaoğlu, *errbudget: choosing the numerical parameters of a black-box solver for a
> target accuracy in one observable*, Computer Physics Communications (Computer Programs
> in Physics), 2026.

It holds the program exactly as the article describes it (version 1.2.1), every analysis
script, the stored solver output of every run, the tables and figures, the pre-registered
held-out tests, and the cluster jobs. It is a record: the program here does not change
except by a new, documented version.

If you want to **use** the method on your own solver rather than check the article, use
the library repository instead: <https://github.com/bioflume/errbudget>. It is the same
engine with a command-line tool, job files and a user guide.

**The name.** The program was called `knobs` until revision 5 of the article (October
2026). Only names changed: the package moved from `src/knobs` to `src/errbudget`, and the
class `Knob` is now `Parameter`. The engine (`_core.py`) is byte-identical, and every table
of the article is reproduced byte for byte (`CHANGELOG.md`). The records fixed by the pre-registration keep the old name:
`PREREGISTRATION.md`, the tags `v1.2.0-frozen` and `v1.2.1`, the frozen copy
`paper/frozen_src/knobs` and the registered J5 driver `hamilton/J5_abc3d/abc_hpc.py`.

## Layout

```
src/errbudget/       the program, version 1.2.1, as the article describes it; do not
                     edit it here (the frozen program of the held-out tests is the tag
                     v1.2.0-frozen, whose hash the replays check)
paper/               everything that reproduces the article on one machine:
                     run_v6_analyses.sh (revision 4), reproduce.sh (revisions 2-3), the
                     analysis scripts, the stored solver output (data/), tables and
                     figures -- start with paper/README.md
hamilton/            everything that runs on Durham's Hamilton cluster (jobs J1-J5,
                     including the pre-registered held-out test H3) -- the runbook is
                     hamilton/README.md
tests/               the behavioural claims of the article as tests (pytest)
examples/            three small examples of the Python API
PREREGISTRATION.md   the held-out tests, fixed before they were run (do not edit)
CHANGELOG.md         what changed between versions, with the comparisons behind it
DESIGN_NOTES.md      why the method is shaped this way, including what failed
THIRD_PARTY_RESULTS.md  the production codes driven through their input files
CITATION.cff         how to cite
```

## Reproducing the article

Three tiers, from cheap to expensive; `paper/README.md` has the details and maps every
number, table and figure to the script that produces it.

| Tier | Command | Needs | Time |
|---|---|---|---|
| A. Replay the stored solver output | `cd paper && bash run_v6_analyses.sh before` | Python 3.9+ with NumPy, SciPy, Matplotlib, SymPy, pytest | 2–4 h on 2 cores |
| B. Rebuild the solver output | `cd paper && ./reproduce.sh`, then `dfg2d2_job/` | tier A plus LAMMPS, GROMACS, ABINIT, OpenFOAM, Quantum ESPRESSO | ~16 CPU-hours + 40 CPU-hours |
| C. Cluster studies J1–J5 | `hamilton/README.md` | a SLURM cluster | 1–3 days |

Quick start for tier A:

```bash
git clone https://github.com/bioflume/errbudget-paper.git && cd errbudget-paper
python -m venv .venv && source .venv/bin/activate
pip install -e ".[paper]"
pytest -q                                   # 52 tests, about 80 s
cd paper && bash run_v6_analyses.sh check   # what the machine has
bash run_v6_analyses.sh before              # the revision-4 analyses
```

Three steps of tier A call a solver (A03 ABINIT, A04 LAMMPS, half of A12 GROMACS); they
fail cleanly without it and the rest goes on. The article's own environment (Python
3.11, NumPy 1.26.4, SciPy 1.13.1, Ubuntu 24.04) is listed in `paper/README.md`; NumPy
2.x and SciPy 1.18 give identical counts.

## The program

```python
from errbudget import Parameter, tune

result = tune(solver,
              [Parameter("N",    [32, 64, 128, 256],        kind="discretisation"),
               Parameter("dt",   [1e-2, 5e-3, 2.5e-3],      kind="discretisation"),
               Parameter("rtol", [1e-4, 1e-6, 1e-8, 1e-10], kind="nested-solver")],
              observable="drag", target=1e-3, start={"N": 64, "dt": 5e-3, "rtol": 1e-10})
print(result)
```

`src/errbudget/api.py` is the public interface (`Parameter`, `tune`, `TuneResult`, `RunFailed`),
`src/errbudget/_core.py` the measurement, the local-order and joint checks and the search,
`src/errbudget/_asymptotic.py` the balance law of the secondary result, and
`src/errbudget/problems/` the test problems with exact solutions. The default configuration is
`errbudget.V121`; `errbudget.V12` is the frozen version of the held-out tests (git tag
`v1.2.0-frozen`, copy in `paper/frozen_src`), `errbudget.V11` the previous version.
`examples/` shows the API on the test problems; the library repository has the user guide.

## Versions and tags

| Tag / version | Meaning |
|---|---|
| `v1.2.0-frozen` | the program frozen before the held-out tests (SHA-256 of the source in `PREREGISTRATION.md`) |
| `v1.2.1` | the released default after the amendment of 11 October 2026; the article's version (package name `knobs`) |
| `main` | this: revision 5 of the article (patch v7) and the rename to `errbudget` |

## Tests

```bash
pip install -e ".[dev]" && pytest
```

The tests encode the article's claims: that recommendations meet their target from a
range of starting points, that unreachable targets are refused, that the diagnostics never
contradict the recommendation, that the measurement recovers an error imposed analytically,
and (`tests/test_spec_v121.py`) the behaviour of version 1.2.1 as Section 4 describes it.

## Citing

See `CITATION.cff`. The archived release with the solver outputs of every run is on Zenodo
(DOI in `CITATION.cff` once minted).

## Licence

MIT.
