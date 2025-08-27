# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Generator
from contextlib import contextmanager
from time import perf_counter_ns
from typing import Callable


@contextmanager
def timer() -> Generator[Callable[[], float], None, None]:
    """Timer context manager."""
    start = perf_counter_ns()
    yield lambda: (perf_counter_ns() - start) / 1e9
