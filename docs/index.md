# errbudget

`errbudget` chooses the numerical parameters of a solver (a mesh, a time step, a cutoff,
the tolerance of a Newton or linear solver inside it) so that one observable reaches a
target accuracy, at close to the cheapest cost that achieves it.

You give it something that runs your solver (a Python function or a shell command), the
numerical parameters with their candidate values, the observable you want, and how
accurate you need it. It gives you a setting, an error estimate, and a statement of which
parameter is limiting the accuracy and which ones you are refining for nothing.

It needs no adjoint, no access to the solver's source, and no analysis of your particular
problem. The number of runs grows with the number of parameters, not with the product of
their candidate values. The tuning itself costs runs, about a hundred to a few hundred in
the article's tests, so it pays back for parameter studies and ensembles but not for a
single run.

## Install

```bash
pip install "errbudget[yaml] @ git+https://github.com/bioflume/errbudget"
```

It needs Python 3.9 or later. The `yaml` extra installs PyYAML, which YAML job files need.

## Quick start

```bash
errbudget init job.yaml    # writes an annotated job file: edit the command, the parameters, the target
errbudget check job.yaml   # validates it and runs your starting setting once
errbudget run job.yaml     # tunes, prints the report, writes errbudget_result.json and .md
```

The [job-file reference](job-files.md) describes every field, and the [recipes](recipes.md)
give job files for OpenFOAM, LAMMPS, GROMACS, Quantum ESPRESSO, ABINIT and MPI codes. For a
solver you can call as a Python function, see the [Python API](api.md). Four
runnable examples, including a command-line one, are in
[`examples/`](https://github.com/bioflume/errbudget/tree/main/examples).

## Contents

```{toctree}
:maxdepth: 1

guide
job-files
recipes
api
changelog
```

## The article

The method and its evaluation are described in an article. Its repository,
<https://github.com/bioflume/errbudget-paper>, has every script, the stored solver output
and the cluster jobs. To cite errbudget, see
[`CITATION.cff`](https://github.com/bioflume/errbudget/blob/main/CITATION.cff).
