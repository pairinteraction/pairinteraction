# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

# ruff: noqa: INP001
from __future__ import annotations

import asyncio
import logging
import os
import shutil
import sys
from functools import partial
from pathlib import Path, PurePath
from subprocess import CalledProcessError
from typing import TYPE_CHECKING, Any

# mypy: disable-error-code="import-untyped"
from sphinx_polyversion.builder import Builder, BuildError
from sphinx_polyversion.driver import DefaultDriver
from sphinx_polyversion.environment import Environment
from sphinx_polyversion.git import Git, closest_tag, refs_by_type
from sphinx_polyversion.json import GLOBAL_ENCODER
from sphinx_polyversion.pyvenv import Pip, VirtualPythonEnvironment
from sphinx_polyversion.sphinx import SphinxBuilder

if TYPE_CHECKING:
    import json
    from types import TracebackType

    from sphinx_polyversion.git import GitRef
    from sphinx_polyversion.json import JSONable

logger = logging.getLogger(__name__)

# Determine repository root directory
DOCS_DIR = Path(__file__).parent
ROOT_DIR = Git.root(DOCS_DIR)


# Define output directory, where the documentation will be built to
OUTPUT_DIR = DOCS_DIR / "_build_polyversion"
shutil.rmtree(OUTPUT_DIR, ignore_errors=True)


# Define all branches and tags to build documentation for (using regex expressions)
LEGACY_VERSION = "v0.9.10"

# for each minor version keep the latest patch version
_git_tags = Git("", tag_regex="v2.*")
_, _all_tags = refs_by_type(asyncio.run(_git_tags.retrieve(ROOT_DIR)))
WANTED_TAGS: list[str] = []
for t in sorted([t.name for t in _all_tags], reverse=True):
    if any(w.startswith(f"v2.{t.split('.')[1]}") for w in WANTED_TAGS):
        continue
    WANTED_TAGS.append(t)
WANTED_TAGS = sorted(WANTED_TAGS)

GIT_OBJ = Git(branch_regex="master", tag_regex="|".join([*WANTED_TAGS, LEGACY_VERSION]), buffer_size=1 * 10**9)
ALL_TARGETS: list[GitRef] = asyncio.run(GIT_OBJ.retrieve(ROOT_DIR))
LATEST = max(t for t in WANTED_TAGS)


# Define factory method for data passed to sphinx,
# important for the version selector (see _templates/versions.html) and the builder (see below)
def data_factory(_driver: DefaultDriver, rev: GitRef, _env: Environment) -> dict[str, Any]:
    version_names = {  # naming convention: "<tag/branch> (<note>)"
        "master": "master (dev)",
        LATEST: f"{LATEST} (stable)",
        "stable": f"{LATEST} (stable)",
        **{tag: tag for tag in WANTED_TAGS[::-1] if tag != LATEST},
        LEGACY_VERSION: f"{LEGACY_VERSION} (legacy)",
    }
    version_tuples = [(name, f"../{key}/index.html") for key, name in version_names.items() if key != LATEST]
    current_version = version_names.get(rev.name, rev.name)
    return {
        "current_rev": rev.name,
        "current_version": current_version,
        "version_tuples": version_tuples,
    }


# Define factory method for data passed to templates, important for the main index.html (see docs/templates/index.html)
def root_data_factory(_driver: DefaultDriver) -> dict[str, str]:
    return {}


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
        apply_patch_if_exists: bool = False,
    ) -> None:
        super().__init__()
        self.cmds = cmds
        self.source = PurePath(source)
        self.logger = logger
        self.encoder = encoder or GLOBAL_ENCODER
        self.pre_cmds = pre_cmds
        self.post_cmds = post_cmds
        self.apply_patch_if_exists = apply_patch_if_exists

    async def build(self, environment: Environment, output_dir: Path, data: JSONable) -> None:  # noqa: C901
        print(f"Building documentation for revision {data['current_rev']}")

        source_dir = str(environment.path.absolute() / self.source)

        def replace(v: str) -> str:
            if "Placeholder" in v:
                v = v.replace("Placeholder.ROOT_DIR", str(environment.path.absolute()))
                v = v.replace("Placeholder.SOURCE_DIR", str(source_dir))
                v = v.replace("Placeholder.OUTPUT_DIR", str(output_dir))
            return v

        patch_cmds: list[tuple[str, ...]] = []
        patch_file = DOCS_DIR / "patches" / f"{data['current_rev']}.patch"
        if self.apply_patch_if_exists and patch_file.exists():
            patch_cmds = [("git", "-C", "Placeholder.ROOT_DIR", "apply", str(patch_file))]

        commands = {
            "patch_cmds": patch_cmds,
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

        for key in ["patch_cmds", "pre_cmds", "cmds", "post_cmds"]:
            cmds = commands[key]
            if cmds is None:
                continue
            for cmd in cmds:
                self.logger.info("Running command (%s): %s", key, " ".join(cmd))
                out, err, rc = await environment.run(*cmd, env=env)
                self.logger.debug("Output:\n %s", out)
                if rc:
                    raise BuildError from CalledProcessError(rc, " ".join(cmd), out, err)


class CustomDriver(DefaultDriver):  # type: ignore [misc]
    def build_failed(self, rev: JSONable, exc_info: tuple[type[BaseException], BaseException, TracebackType]) -> None:
        """Override to exit with error code 1 on build failure."""
        super().build_failed(rev, exc_info)
        sys.exit(1)


# Mapping of versions/revisions to builders and environments, which is used for building the documentation
# Versions/Revisions not listed here will use the entry, which revision is the closest ancestor of the wanted revision
# IMPORTANT: The revisions must be put in in the correct order (starting with the oldest)
BUILDER_MAPPING = {
    LEGACY_VERSION: CustomCommandBuilder(
        Path("doc/sphinx/"),
        cmds=[("sphinx-build", "--color", "-av", "--keep-going", "Placeholder.SOURCE_DIR", "Placeholder.OUTPUT_DIR")],
        apply_patch_if_exists=True,
    ),
    # apparently we need a separate builder instance for each tag, otherwise sphinx-polyversion will fail
    **{
        tag: CustomCommandBuilder(
            Path("docs/"),
            cmds=[
                ("sphinx-build", "--color", "-avW", "--keep-going", "Placeholder.SOURCE_DIR", "Placeholder.OUTPUT_DIR")
            ],
            apply_patch_if_exists=True,
        )
        for tag in WANTED_TAGS
    },
    "master": SphinxBuilder(Path("docs/"), args=["-a", "-v", "-W", "--keep-going"]),
}

ENVIRONMENT_MAPPING = {
    LEGACY_VERSION: Pip.factory(
        venv=ROOT_DIR / "tmp_venv",
        args=[
            f"pairinteraction=={LEGACY_VERSION[1:]}",
            "sphinx<9",
            "nbsphinx",
            "sphinxcontrib-mermaid<1",
            "sphinx_rtd_theme",
            "sphinx_polyversion",
        ],
    ),
    **{
        tag: Pip.factory(
            venv=ROOT_DIR / "tmp_venv",
            args=[f"pairinteraction[docs]=={tag[1:]}"],
        )
        for tag in WANTED_TAGS
    },
    "master": VirtualPythonEnvironment.factory(venv=ROOT_DIR / ".venv"),
}

# Create the actual driver instance and run it to build all wanted documentations
driver = CustomDriver(
    ROOT_DIR,
    output_dir=OUTPUT_DIR,
    vcs=GIT_OBJ,
    builder=BUILDER_MAPPING,
    env=ENVIRONMENT_MAPPING,
    selector=partial(closest_tag, ROOT_DIR),
    template_dir=DOCS_DIR / "templates",
    static_dir=DOCS_DIR / "static",
    root_data_factory=root_data_factory,
    data_factory=data_factory,
)
driver.run(sequential=True)

# rename latest "v2.x.x/" folder to "stable/"
stable_dir = OUTPUT_DIR / LATEST
if stable_dir.exists():
    shutil.move(stable_dir, OUTPUT_DIR / "stable")
