<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/errbudget-logo-dark.svg">
    <img alt="errbudget" src="docs/assets/errbudget-logo.svg" width="360">
  </picture>
</p>

# Reproducing the errbudget article

This branch reproduces the article "errbudget: choosing the discretization and the solver
tolerances of a black-box code for a target accuracy in one observable" (Computer Physics
Communications) and its supplement. The program for users is on the branch `main`.

The article presents only the released program, errbudget 1.2.1, with its default
configuration `errbudget.V121`, and every result was replayed with it. The configurations
`V12` and `V11` and the tag `v1.2.0-frozen` exist only to replay the pre-registered searches
of H1–H3, which ran the program as frozen at registration. The program was earlier named
`knobs`; the frozen source keeps that name.

## What is here

```
src/errbudget/   the package as it produced the results: the engine (_core.py), the Python
                 interface (api.py: tune, Parameter) and the test problems
paper/           the replay and rebuild scripts, the stored runs (data/), the tables and the
                 printed output (figures/), the frozen program (frozen_src/) and the build of
                 the coarser 2D-2 study (dfg2d2_job/)
hamilton/        the Hamilton jobs J1-J5 and their summarised results (results/)
studies/         pre_comparison/ (Section S-4.6), reuse2d1/ (S-6), airfoil/ (S-7)
tests/           the tests, among them the 14 of test_spec_v121.py (Section S-1.11)
examples/        short examples of the Python interface
```

`PREREGISTRATION.md` fixed the pre-registered problems H1–H3. `DATA.md` lists the raw data
kept outside git, with their SHA-256. `CHANGELOG.md` records how the repository changed.

## Reproduce

```bash
git checkout paper
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e . -r paper/requirements-article.txt
bash paper/replay_all.sh
```

`paper/README.md` maps every number, table and figure of the article and the supplement to
the script that produces it, and says which results need a solver (`paper/reproduce.sh`) or
a cluster (`hamilton/`, `studies/`). `paper/REPRODUCTION_CHECK.md` records what was
reproduced.

## The command-line tool

The commands `errbudget init`, `check` and `run` of Section 9 belong to the package on the
branch `main` (`src/errbudget/cli.py`, installed as the console script `errbudget`). The
package here is the record of the code that produced the results; it has the engine and
the Python interface but not the command-line tool.

## Licence and citation

MIT (`LICENSE`). To cite the article or the software, see `CITATION.cff`.
