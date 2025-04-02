# SPDX-FileCopyrightText: 202 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import os
from typing import Union

from pairinteraction._backend import run_unit_tests


def run_module_tests(download_missing: bool = False, database_dir: Union[str, os.PathLike[str]] = "") -> int:
    use_cache = True
    return run_unit_tests(download_missing, use_cache, database_dir)
