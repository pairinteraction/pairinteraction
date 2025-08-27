# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: INP001

import json
import logging
import os
import shutil
from functools import partial
from pathlib import Path, PurePath
from subprocess import CalledProcessError
from typing import Any, Union

# mypy: disable-error-code="import-untyped"
from sphinx_polyversion.builder import Builder, BuildError
from sphinx_polyversion.driver import DefaultDriver
from sphinx_polyversion.environment import Environment
from sphinx_polyversion.git import Git, GitRef, closest_tag, refs_by_type
from sphinx_polyversion.json import GLOBAL_ENCODER, JSONable
from sphinx_polyversion.pyvenv import VirtualPythonEnvironment
from sphinx_polyversion.sphinx import SphinxBuilder

import pairinteraction

logger = logging.getLogger(__name__)

# Determine repository root directory
DOCS_DIR = Path(__file__).parent
ROOT_DIR = Git.root(DOCS_DIR)


# Define output directory, where the documentation will be built to
OUTPUT_DIR = DOCS_DIR / "_build_polyversion"
shutil.rmtree(OUTPUT_DIR, ignore_errors=True)


# Define all branches and tags to build documentation for (using regex expressions)
LATEST_BRANCH = "master"
LEGACY_VERSION = "v0.9.10"


# Define factory method for data passed to sphinx, important for the version selector (see _templates/versions.html)
def data_factory(_driver: DefaultDriver, rev: GitRef, _env: Environment) -> dict[str, Any]:
    version_names = {
        LATEST_BRANCH: f"latest (v{pairinteraction.__version__})",
        LEGACY_VERSION: f"legacy ({LEGACY_VERSION})",
    }
    base_url = "../"
    version_tuples = [
        [version_names[LATEST_BRANCH], base_url + "index.html"],
        [version_names[LEGACY_VERSION], base_url + f"{LEGACY_VERSION}/index.html"],
    ]

    current_version = version_names.get(rev.name, rev.name)
    return {
        "current_version": current_version,
        "version_tuples": version_tuples,
    }


# Define factory method for data passed to templates, important for the main index.html (see docs/templates/index.html)
def root_data_factory(driver: DefaultDriver) -> dict[str, Union[list[GitRef], GitRef, None]]:
    revisions: list[GitRef] = driver.builds
    branches, _tags = refs_by_type(revisions)
    latest = next((b for b in branches if b.name == LATEST_BRANCH), None)
    return {"revisions": revisions, "latest": latest}


# Define a custom documentation Builder
# (mainly allow for replacing placeholders in strings and running multiple commands)
class CustomCommandBuilder(Builder[Environment, None]):  # type: ignore [misc]
    def __init__(
        self,
        source: str | PurePath,
        cmds: list[tuple[str, ...]],
        encoder: json.JSONEncoder | None = None,
        pre_cmds: list[tuple[str, ...]] | None = None,
        post_cmds: list[tuple[str, ...]] | None = None,
    ) -> None:
        super().__init__()
        self.cmds = cmds
        self.source = PurePath(source)
        self.logger = logger
        self.encoder = encoder or GLOBAL_ENCODER
        self.pre_cmds = pre_cmds
        self.post_cmds = post_cmds

    async def build(self, environment: Environment, output_dir: Path, data: JSONable) -> None:
        self.logger.info("Building...")

        source_dir = str(environment.path.absolute() / self.source)

        def replace(v: str) -> str:
            if "Placeholder" in v:
                v = v.replace("Placeholder.OUTPUT_DIR", str(output_dir))
                v = v.replace("Placeholder.SOURCE_DIR", str(source_dir))
            return v

        commands = {
            "pre_cmds": self.pre_cmds,
            "cmds": self.cmds,
            "post_cmds": self.post_cmds,
        }

        for cmds in commands.values():
            if cmds is None:
                continue
            for i, cmd in enumerate(cmds):
                cmds[i] = tuple(map(replace, cmd))

        env = os.environ.copy()
        env["POLYVERSION_DATA"] = self.encoder.encode(data)

        # create output directory
        output_dir.mkdir(exist_ok=True, parents=True)

        for key in ["pre_cmds", "cmds", "post_cmds"]:
            cmds = commands[key]
            if cmds is None:
                continue
            for cmd in cmds:
                self.logger.info("Running command (%s): %s", key, " ".join(cmd))
                out, err, rc = await environment.run(*cmd, env=env)
                self.logger.debug("Output:\n %s", out)
                if rc:
                    raise BuildError from CalledProcessError(rc, " ".join(cmd), out, err)


# Mapping of versions/revisions to builders and environments, which is used for building the documentation
# Versions/Revisions not listed here will use the entry, which revision is the closest ancestor of the wanted revision
# IMPORTANT: The revisions must be put in in the correct order (starting with the oldest)
BUILDER_MAPPING = {
    LEGACY_VERSION: CustomCommandBuilder(
        Path("doc/sphinx/"),
        pre_cmds=[("git", "-C", "Placeholder.SOURCE_DIR/../..", "apply", f"{DOCS_DIR}/patches/{LEGACY_VERSION}.patch")],
        cmds=[
            ("sphinx-build", "--color", "-a", "-v", "--keep-going", "Placeholder.SOURCE_DIR", "Placeholder.OUTPUT_DIR")
        ],
    ),
    "v2.0.0": SphinxBuilder(Path("docs/"), args=["-a", "-v", "-W", "--keep-going"]),
}

ENVIRONMENT_MAPPING = {
    LEGACY_VERSION: VirtualPythonEnvironment.factory(venv=ROOT_DIR / "venv_pairinteraction_v0.9.10"),
    "v2.0.0": VirtualPythonEnvironment.factory(venv=ROOT_DIR / ".venv"),
}


# Create the actual driver instance and run it to build all wanted documentations
DefaultDriver(
    ROOT_DIR,
    output_dir=OUTPUT_DIR,
    vcs=Git(
        branch_regex=LATEST_BRANCH,
        tag_regex=LEGACY_VERSION,
        buffer_size=1 * 10**9,  # 1 GB
    ),
    builder=BUILDER_MAPPING,
    env=ENVIRONMENT_MAPPING,
    selector=partial(closest_tag, ROOT_DIR),
    template_dir=DOCS_DIR / "templates",
    static_dir=DOCS_DIR / "static",
    root_data_factory=root_data_factory,
    data_factory=data_factory,
).run(sequential=True)
