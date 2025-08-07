# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: INP001

import shutil
from functools import partial
from pathlib import Path
from typing import Any, Union

# mypy: disable-error-code="import-untyped"
from sphinx_polyversion.driver import DefaultDriver
from sphinx_polyversion.environment import Environment
from sphinx_polyversion.git import Git, GitRef, closest_tag, refs_by_type
from sphinx_polyversion.json import JSONable
from sphinx_polyversion.pyvenv import VirtualPythonEnvironment
from sphinx_polyversion.sphinx import CommandBuilder, Placeholder, SphinxBuilder

import pairinteraction

# Determine repository root directory
DOCS_DIR = Path(__file__).parent
ROOT_DIR = Git.root(DOCS_DIR)


# Define output directory, where the documentation will be built to
OUTPUT_DIR = DOCS_DIR / "_build_polyversion"
shutil.rmtree(OUTPUT_DIR, ignore_errors=True)


# Define all branches and tags to build documentation for (using regex expressions)
BRANCH_REGEX = r"master"
TAG_REGEX = r"v0.9.10"


# Define factory method for data passed to sphinx, important for the version selector (see _templates/versions.html)
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


# Define factory method for data passed to templates, important for the main index.html (see docs/templates/index.html)
def root_data(driver: DefaultDriver) -> dict[str, Union[list[GitRef], GitRef, None]]:
    revisions: list[GitRef] = driver.builds
    branches, _tags = refs_by_type(revisions)
    latest = next((b for b in branches if b.name == "master"), None)
    return {"revisions": revisions, "latest": latest}


# Define a custom documentation Builder (mainly allow for replacing placeholders in strings)
class UpdatedCommandBuilder(CommandBuilder):  # type: ignore [misc]
    pre_cmd: tuple[Any, ...]
    cmd: tuple[Any, ...]
    post_cmd: tuple[Any, ...]

    async def build(self, environment: Environment, output_dir: Path, data: JSONable) -> None:
        source_dir = str(environment.path.absolute() / self.source)

        def replace(v: Any) -> Any:
            if isinstance(v, str) and "Placeholder" in v:
                v = v.replace("Placeholder.OUTPUT_DIR", str(output_dir))
                v = v.replace("Placeholder.SOURCE_DIR", str(source_dir))
            return v

        if self.pre_cmd:
            self.pre_cmd = tuple(map(replace, self.pre_cmd))
        if self.cmd:
            self.cmd = tuple(map(replace, self.cmd))
        if self.post_cmd:
            self.post_cmd = tuple(map(replace, self.post_cmd))
        await super().build(environment, output_dir, data)


# Mapping of versions/revisions to builders and environments, which is used for building the documentation
# Versions/Revisions not listed here will use the entry, which revision is the closest ancestor of the wanted revision
# IMPORTANT: The revisions must be put in in the correct order (starting with the oldest)
BUILDER_MAPPING = {
    "v0.9.10": UpdatedCommandBuilder(
        Path("doc/sphinx/"),
        pre_cmd=["cp", "Placeholder.SOURCE_DIR/conf.py.in", "Placeholder.SOURCE_DIR/conf.py"],
        cmd=["sphinx-build", "--color", "-a", "-v", "--keep-going", Placeholder.SOURCE_DIR, Placeholder.OUTPUT_DIR],
    ),
    "v2.0.0": SphinxBuilder(Path("docs/"), args=["-a", "-v", "-W", "--keep-going"]),
}

ENVIRONMENT_MAPPING = {
    "v0.9.10": VirtualPythonEnvironment.factory(venv=ROOT_DIR / "venv_pairinteraction_v0.9.10"),
    "v2.0.0": VirtualPythonEnvironment.factory(venv=ROOT_DIR / ".venv"),
}


# Create the actual driver instance and run it to build all wanted documentations
DefaultDriver(
    ROOT_DIR,
    output_dir=OUTPUT_DIR,
    vcs=Git(
        branch_regex=BRANCH_REGEX,
        tag_regex=TAG_REGEX,
        buffer_size=1 * 10**9,  # 1 GB
    ),
    builder=BUILDER_MAPPING,
    env=ENVIRONMENT_MAPPING,
    selector=partial(closest_tag, ROOT_DIR),
    template_dir=DOCS_DIR / "templates",
    static_dir=DOCS_DIR / "static",
    root_data_factory=root_data,
    data_factory=data,
).run(sequential=True)
