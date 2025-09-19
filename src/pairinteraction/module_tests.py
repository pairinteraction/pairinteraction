# SPDX-FileCopyrightText: 202 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING

from pairinteraction._backend import run_unit_tests

if TYPE_CHECKING:
    import os


def run_module_tests(download_missing: bool = False, database_dir: str | os.PathLike[str] = "") -> int:
    use_cache = True
    return run_unit_tests(download_missing, use_cache, database_dir)
