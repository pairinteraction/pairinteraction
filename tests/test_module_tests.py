"""
Execute the module tests.
"""

from pairinteraction import run_module_tests


def test_module_tests(download_missing: bool, database_dir: str) -> None:
    """Execute the module tests."""
    run_module_tests(download_missing, database_dir)
