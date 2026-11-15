## Code style

- Write what an experienced research software engineer would keep after review: every line has a job.
- The engine `src/errbudget/_core.py` stays byte-identical to the one that produced the article's results (sha256 55f0e994...) until engine 1.2.2; its planned clean-up is in `CHANGELOG.md` under "Unreleased".
- The package is `src/errbudget/`: `api.py` (`tune`, `Parameter`), `jobfile.py` and `runner.py` (job files and runs of a command), `cli.py` (`errbudget init/check/run`), `report.py`, `problems/` (the test problems). Examples in `examples/` import the package; they never edit `sys.path` or copy code.
- Validate input where it enters (`Parameter`, job files, the command line) with a clear message; not on internal calls.
- No `try/except` unless the handler does something specific (report a failed run, retry a parse). Never swallow an error or fall back silently.
- Comments explain why (equation, rule of the article, numerical reason), never what. No banners, emoji or commented-out code.
- Numpydoc docstrings on the public functions (`tune`, `Parameter`, `run_protocol`, `load_job`, `run_job`), with what each returns; one line or none on private helpers.
- No logging configuration or prints inside the package, except the command-line tool and `verbose=True`.
- A new option of the program is a key of the configuration (`V121`) or of the job file, documented in `docs/`, with a test in `tests/`.
- Run `pytest tests -q` before every commit; the tests of `tests/test_spec_v121.py` pin the rules of the article.
