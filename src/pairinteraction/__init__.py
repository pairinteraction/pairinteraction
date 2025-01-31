import os

from pairinteraction import backend
from pairinteraction.backend._backend import VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH
from pairinteraction.module_tests import run_module_tests
from pairinteraction.setup_test_mode import setup_test_mode
from pairinteraction.units import ureg

__all__ = [
    "backend",
    "run_module_tests",
    "ureg",
]

__version__ = f"{VERSION_MAJOR}.{VERSION_MINOR}.{VERSION_PATCH}"


if os.getenv("PAIRINTERACTION_TEST_MODE", "0") == "1":
    setup_test_mode()
