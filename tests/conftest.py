from pathlib import Path

import pytest
from pint import UnitRegistry


def pytest_addoption(parser) -> None:
    parser.addoption("--generate-reference", action="store_true", default=False, help="Generate reference data")
    parser.addoption(
        "--database-dir",
        action="store",
        default=str(Path(__file__).parent.parent / "data/database"),
        help="Path to the database directory",
    )


@pytest.fixture(scope="session")
def generate_reference(pytestconfig) -> bool:
    return pytestconfig.getoption("--generate-reference")


@pytest.fixture(scope="session")
def database_dir(pytestconfig) -> Path:
    return Path(pytestconfig.getoption("--database-dir"))


@pytest.fixture(scope="session")
def ureg() -> UnitRegistry:
    """Create and return a UnitRegistry with atomic units."""
    return UnitRegistry(system="atomic")
