# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: INP001

from pathlib import Path
from typing import Any, Union

from sphinx_polyversion.api import apply_overrides
from sphinx_polyversion.driver import DefaultDriver
from sphinx_polyversion.environment import Environment
from sphinx_polyversion.git import Git, GitRef, refs_by_type
from sphinx_polyversion.pyvenv import Pip
from sphinx_polyversion.sphinx import SphinxBuilder

import pairinteraction

#: Regex matching the branches to build docs for
BRANCH_REGEX = r"master|docs-polyversion"  # FIXME docs-polyversion just for testing

#: Regex matching the tags to build docs for
TAG_REGEX = r"v2.*"
# TAG_REGEX = r"v0.9.10"  # v0.9.10 has different docs structure

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


# Data passed to sphinx, important for the version selector
def data(_driver: DefaultDriver, rev: GitRef, _env: Environment) -> dict[str, Any]:
    base_url = "https://www.pairinteraction.org/pairinteraction/sphinx/html/"
    versions_dict = [
        [f"latest ({pairinteraction.__version__})", base_url],
        ["old pairinteraction (v0.9.10)", base_url + "v0.9.10/"],
    ]

    current_version = rev.name
    if current_version == "master":
        current_version = f"latest ({pairinteraction.__version__})"
    return {
        "current_version": current_version,
        "versions_dict": versions_dict,
    }


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
    ),
    builder=SphinxBuilder(src, args=SPHINX_ARGS),
    env=Pip.factory(venv="../.venv", args=PIP_ARGS),
    template_dir=root / src / "templates",
    static_dir=root / src / "static",
    root_data_factory=root_data,
    data_factory=data,
).run(sequential=False)
