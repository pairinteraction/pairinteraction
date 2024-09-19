from pathlib import Path
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
        default=str(Path(__file__).parent.parent / "data/database"),
        help="Path to the database directory",
    )
    parser.addoption("--download-missing", action="store_true", default=False, help="Download missing database files")


@pytest.fixture(scope="session")
def generate_reference(pytestconfig: "Config") -> bool:
    return pytestconfig.getoption("--generate-reference")


@pytest.fixture(scope="session")
def database_dir(pytestconfig: "Config") -> str:
    return pytestconfig.getoption("--database-dir")


@pytest.fixture(scope="session")
def download_missing(pytestconfig: "Config") -> bool:
    return pytestconfig.getoption("--download-missing")


@pytest.fixture(scope="session")
def ureg() -> UnitRegistry:
    """Create and return a UnitRegistry with atomic units."""
    return UnitRegistry(system="atomic")
