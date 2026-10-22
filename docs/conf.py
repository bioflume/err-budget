# Configuration of the errbudget documentation (Sphinx, MyST-Parser, Furo).
#
# From the repository root:
#     pip install -r docs/requirements.txt
#     sphinx-build -b html docs docs/_build/html

import re
from pathlib import Path

# -- Project -------------------------------------------------------------------

project = "errbudget"
author = "the errbudget authors"
copyright = "2026, the errbudget authors"


def _version_from_pyproject():
    """Return the version in the [project] table of ../pyproject.toml.

    It is read with a regular expression, not with tomllib (Python 3.11+) or by
    importing errbudget, so that the documentation builds on Python 3.10 without
    the package installed.
    """
    text = (Path(__file__).resolve().parent.parent / "pyproject.toml").read_text(encoding="utf-8")
    table = re.search(r"^\[project\][ \t]*$(.*?)(?=^\[|\Z)", text, re.MULTILINE | re.DOTALL)
    match = table and re.search(r"""^version\s*=\s*["']([^"']+)["']""", table.group(1), re.MULTILINE)
    if not match:
        raise RuntimeError('docs/conf.py: no version = "..." in the [project] table of pyproject.toml')
    return match.group(1)


release = _version_from_pyproject()
version = release

# -- General -------------------------------------------------------------------

extensions = ["myst_parser"]

# Every .md file in docs/ is a page, except a README.md (GitHub shows one when the
# folder is browsed). docs/assets/ is skipped already, as html_static_path.
exclude_patterns = ["_build", "README.md", "**/README.md"]

# Code blocks without a language hold solver input or report text, not Python.
highlight_language = "none"

# Anchors for headings down to ###, so that links such as #exit-codes or
# job-files.md#exit-codes work.
myst_heading_anchors = 3

# -- HTML output ---------------------------------------------------------------

html_theme = "furo"
html_title = "errbudget"
html_static_path = ["assets"]
html_favicon = "assets/errbudget-tile-32.png"
html_theme_options = {
    # File names relative to html_static_path, i.e. in docs/assets/.
    "light_logo": "errbudget-logo.svg",
    "dark_logo": "errbudget-logo-dark.svg",
    # The logo carries the name.
    "sidebar_hide_name": True,
}
