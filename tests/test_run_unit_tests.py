# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Execute the module tests."""

import pairinteraction as pi
from pairinteraction import run_unit_tests


def test_module_tests() -> None:
    """Execute the module tests."""
    database = pi.Database.get_global_database()
    assert run_unit_tests(database.download_missing, database_dir=database.database_dir) == 0
