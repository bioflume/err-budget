# Contributing to errbudget

Thank you for considering a contribution. errbudget has two parts, and they are held to
different rules. The engine, `src/errbudget/_core.py`, is the code that the article in
Computer Physics Communications describes, and it does not change without a new engine
version. Everything around it (the command-line interface, job files, the runner, the
report, the documentation and the examples) changes in the ordinary way.

## Setting up

Fork the repository on GitHub, clone your fork, and install the package in editable mode
with the development extras, preferably in a virtual environment:

```bash
git clone https://github.com/YOUR-USERNAME/errbudget
cd errbudget
pip install -e ".[dev,yaml]"
```

You need Python 3.9 or later. The `dev` extra installs pytest, PyYAML and SymPy (some of
the test problems in `errbudget.problems` are manufactured solutions built with SymPy),
and `yaml` is the extra that users install to read YAML job files.

## Running the tests

```bash
pytest                          # the whole suite, a minute or two
pytest tests/test_jobfile.py    # one file
```

The tests encode the article's claims, not only the health of the code.
`tests/test_tune.py` checks recommended settings against exact solutions, and
`tests/test_spec_v121.py` pins the behaviour of the engine that Section 4 of the article
specifies. One test is skipped unless `paper/data/burgers_tensor.pkl`, stored solver
output from the article, is present. A fresh clone does not have it, so the skip is
expected.

## Running the command-line example

`examples/04_command_line` tunes a small solver that runs as a command, through a job
file, as a user would. CI runs the same two commands:

```bash
cd examples/04_command_line
errbudget check job.yaml   # validates the job file and runs the starting setting once
errbudget run job.yaml     # tunes in a few seconds and prints the report
```

The run writes `errbudget_runs/`, `errbudget_result.json` and `errbudget_result.md` into
that directory; git ignores them.

## Reporting a bug or asking for a feature

Open an issue at <https://github.com/bioflume/errbudget/issues/new/choose> and choose
**Bug report** or **Feature request**. The bug form asks for what you ran (the job file or
the code), what happened, what you expected, the output of `errbudget version`, and your
Python version and operating system. A small case that shows the problem is the most
useful thing you can send. If your solver cannot be shared, say what it is and how
errbudget calls it.

A result that declines to certify a setting is not necessarily a bug. errbudget declines
when it cannot bound a parameter's error, and the report says which parameter; extending
that parameter's list of values is often the answer. The
[documentation](https://bioflume.github.io/errbudget/) explains when and why it declines.

Do not report a security problem in a public issue; see the
[security policy](.github/SECURITY.md).

## The engine rule

The engine is version 1.2.1 (`errbudget.ENGINE_VERSION`), the version the article
describes, and `_core.py` is kept byte for byte as that version. A change to it changes
results that the article reports, and this holds for bug fixes too. Such a change
therefore needs

1. **an issue first**, describing the change and why it is needed, and agreement in that
   issue before a pull request is opened;
2. **a new `ENGINE_VERSION`** in `src/errbudget/__init__.py`;
3. **an entry in `CHANGELOG.md`** that says how the results change, not only how the code
   changes;
4. **`tests/test_spec_v121.py` updated or extended**, so that its tests state the new
   behaviour.

Where possible, put new behaviour behind an option, so that a configuration still
reproduces the article, as `errbudget.V12` and `errbudget.V11` reproduce the earlier
versions. Cosmetic changes to `_core.py` (formatting, comments, renaming, including the
engine's own uses of "knob") are not accepted, since they break the byte-for-byte match
without changing anything for users.

A change outside `_core.py` needs none of this. A small fix can go straight to a pull
request; for a larger change, an issue first saves work on both sides.

## Pull requests

Before you open one, check that

- `pytest` passes;
- `CHANGELOG.md` has an entry for any change a user would notice, at the top under
  `## Unreleased` (add the heading if it is missing);
- the engine, `src/errbudget/_core.py`, is unchanged, or the engine rule above is
  followed;
- the documentation (`README.md`, `docs/`, docstrings) describes the new behaviour.

The pull-request template repeats this list. CI runs the tests on Python 3.9, 3.11, 3.13
and 3.14, and runs the command-line example. No formatter is configured; follow the style
of the code around your change.

## Code of conduct

Everyone who takes part in the project is expected to follow its
[code of conduct](.github/CODE_OF_CONDUCT.md).

## Licence

errbudget is released under the [MIT licence](LICENSE). By contributing, you agree that
your contribution is released under the same licence.
