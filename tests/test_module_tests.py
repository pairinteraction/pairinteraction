"""Execute the module tests."""

import pairinteraction.backend.float as pi
from pairinteraction import run_module_tests


def test_module_tests() -> None:
    """Execute the module tests."""
    database = pi.Database.get_global_database()
    assert run_module_tests(database.download_missing, database.database_dir) == 0
