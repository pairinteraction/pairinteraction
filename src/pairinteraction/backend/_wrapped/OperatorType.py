from typing import Literal

from pairinteraction.backend import _backend

OperatorType = Literal[
    "ZERO",
    "ENERGY",
    "ELECTRIC_DIPOLE",
    "ELECTRIC_QUADRUPOLE",
    "ELECTRIC_QUADRUPOLE_ZERO",
    "ELECTRIC_OCTUPOLE",
    "MAGNETIC_DIPOLE",
    "ARBITRARY",
]

CPPOperatorType = _backend.OperatorType
