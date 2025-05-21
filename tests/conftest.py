# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import multiprocessing
import os
from typing import TYPE_CHECKING, Optional

import pytest
from pint import UnitRegistry

if TYPE_CHECKING:
    from _pytest.config import Config
    from _pytest.config.argparsing import Parser
    from pairinteraction_gui.app import Application


def pytest_addoption(parser: "Parser") -> None:
    parser.addoption("--generate-reference", action="store_true", default=False, help="Generate reference data")
    parser.addoption(
        "--database-dir",
        action="store",
        default=None,
        help="Path to the database directory",
    )
    parser.addoption("--download-missing", action="store_true", default=False, help="Download missing database files")


@pytest.fixture(scope="session")
def generate_reference(pytestconfig: "Config") -> bool:
    return pytestconfig.getoption("--generate-reference")  # type: ignore [no-any-return]


@pytest.fixture(scope="session")
def ureg() -> UnitRegistry:
    """Create and return a UnitRegistry with atomic units."""
    return UnitRegistry(system="atomic")


def pytest_sessionstart(session: pytest.Session) -> None:
    """Initialize everything before the tests are run."""
    download_missing: bool = session.config.getoption("--download-missing")
    database_dir: Optional[str] = session.config.getoption("--database-dir")

    # Disable the test mode of pairinteraction that would call _setup_test_mode
    # automatically. This would be necessary for testing the jupyter notebooks
    # but we want to call _setup_test_mode manually.
    test_mode = os.environ.get("PAIRINTERACTION_TEST_MODE", "1")
    os.environ["PAIRINTERACTION_TEST_MODE"] = "0"

    # Call _setup_test_mode manually with the given options
    from pairinteraction import _setup_test_mode

    _setup_test_mode(download_missing, database_dir)

    # Set the test mode environment variables, so they can be used by subprocesses
    os.environ["PAIRINTERACTION_TEST_MODE"] = test_mode
    if database_dir is not None:
        os.environ["PAIRINTERACTION_TEST_DOWNLOAD_MISSING"] = str(int(download_missing))
        os.environ["PAIRINTERACTION_TEST_DATABASE_DIR"] = database_dir

    # For pairinteraction_gui set the multiprocessing start method to "spawn" (see also pairinteraction_gui/__init__.py)
    multiprocessing.set_start_method("spawn")


@pytest.fixture(scope="session")
def qapp_cls() -> type["Application"]:
    """Let the qapp and qtbot fixtures use our custom Application class."""
    from pairinteraction_gui.app import Application

    return Application
