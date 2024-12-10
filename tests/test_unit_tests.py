"""
Execute the unit tests.
"""

import pairinteraction.backend._backend as pi  # TODO access as pairinteraction.run_unit_tests


def test_unit_tests(database_dir: str, download_missing: bool) -> None:
    """Execute the unit tests."""
    pi.run_unit_tests(download_missing, True, database_dir)
