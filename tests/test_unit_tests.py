"""
Execute the unit tests.
"""

import pairinteraction.backend._backend as pi  # TODO access as pairinteraction.backend.double.test


def test_unit_tests(database_dir: str, download_missing: bool) -> None:
    """Execute the unit tests."""
    pi.test(download_missing, True, database_dir)
