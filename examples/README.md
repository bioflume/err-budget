# Examples

| | what it shows | run |
|---|---|---|
| `01_quickstart.py` | `tune` on a solver with an exact solution (1D Burgers), and a check of the advice against the exact value | `python examples/01_quickstart.py` |
| `02_bring_your_own_solver.py` | an ordinary function with no knowledge of knobs: a stiff ODE with the three kinds of parameter in miniature; 15x cheaper than the cautious start | `python examples/02_bring_your_own_solver.py` |
| `03_declaration_falsified.py` | what happens when an influential parameter is declared `cost-only`: the preliminary sweep (`screen=True`) contradicts the declaration, and the recommendation is still safe | `python examples/03_declaration_falsified.py` |
| `04_command_line/` | the same stiff ODE as a command-line program tuned through a job file, no Python written | `cd examples/04_command_line && knobs check job.yaml && knobs run job.yaml` |

The first three need `pip install -e .` from the repository root; the fourth also needs
PyYAML (`pip install -e ".[yaml]"`). Job files for OpenFOAM, LAMMPS, GROMACS and
plane-wave DFT codes are in `docs/recipes.md`.
