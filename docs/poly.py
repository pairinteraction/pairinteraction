# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: INP001

from pathlib import Path
from typing import Union

from sphinx_polyversion.api import apply_overrides
from sphinx_polyversion.driver import DefaultDriver
from sphinx_polyversion.git import Git, GitRef, file_predicate, refs_by_type
from sphinx_polyversion.pyvenv import Pip
from sphinx_polyversion.sphinx import SphinxBuilder

#: Regex matching the branches to build docs for
BRANCH_REGEX = r"docs-polyversion"

#: Regex matching the tags to build docs for
TAG_REGEX = r""
TAG_REGEX = r"v2\.2\.1|v0\.9\.10"

#: Output dir relative to project root
OUTPUT_DIR = "_build"

#: Source directory
SOURCE_DIR = "docs/"

#: Arguments to pass to `pip install`
PIP_ARGS = ["sphinx"]

#: Arguments to pass to `sphinx-build`
SPHINX_ARGS = ["-a", "-v"]
SPHINX_ARGS += ["-W", "--keep-going"]

# Load overrides read from commandline to global scope
apply_overrides(globals())
# Determine repository root directory
root = Git.root(Path(__file__).parent)

# Setup driver and run it
src = Path(SOURCE_DIR)


# Data passed to templates, important for the root index.html (see docs/templates/index.html)
def root_data(driver: DefaultDriver) -> dict[str, Union[list[GitRef], GitRef, None]]:
    revisions: list[GitRef] = driver.builds
    branches, _tags = refs_by_type(revisions)
    latest = next((b for b in branches if b.name == "master"), None)
    return {"revisions": revisions, "latest": latest}


DefaultDriver(
    root,
    OUTPUT_DIR,
    vcs=Git(
        branch_regex=BRANCH_REGEX,
        tag_regex=TAG_REGEX,
        buffer_size=1 * 10**9,  # 1 GB
        predicate=file_predicate([src]),  # exclude refs without source dir
    ),
    builder=SphinxBuilder(src, args=SPHINX_ARGS),
    env=Pip.factory(venv="../.venv", args=PIP_ARGS),
    template_dir=root / src / "templates",
    static_dir=root / src / "static",
    root_data_factory=root_data,
).run(sequential=False)
