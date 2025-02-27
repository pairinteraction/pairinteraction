import os
from typing import TYPE_CHECKING

import pytest
from pint import UnitRegistry

if TYPE_CHECKING:
    from _pytest.config import Config
    from _pytest.config.argparsing import Parser


def pytest_addoption(parser: "Parser") -> None:
    parser.addoption("--generate-reference", action="store_true", default=False, help="Generate reference data")
    parser.addoption(
        "--database-dir",
        action="store",
        default="",
        help="Path to the database directory",
    )
    parser.addoption("--download-missing", action="store_true", default=False, help="Download missing database files")


@pytest.fixture(scope="session")
def generate_reference(pytestconfig: "Config") -> bool:
    return pytestconfig.getoption("--generate-reference")


@pytest.fixture(scope="session")
def ureg() -> UnitRegistry:
    """Create and return a UnitRegistry with atomic units."""
    return UnitRegistry(system="atomic")


def pytest_sessionstart(session: pytest.Session) -> None:
    """Initialize everything before the tests are run."""
    download_missing = session.config.getoption("--download-missing")
    database_dir = session.config.getoption("--database-dir")

    # Make it possible to overwrite the database in pytest test mode
    test_mode = os.getenv("PAIRINTERACTION_TEST_MODE", "1")
    os.environ["PAIRINTERACTION_TEST_MODE"] = "0"
    from pairinteraction import _setup_test_mode

    _setup_test_mode(download_missing, database_dir)
    os.environ["PAIRINTERACTION_TEST_MODE"] = test_mode
