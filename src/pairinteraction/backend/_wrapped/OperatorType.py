from typing import Literal, get_args

from pairinteraction.backend import _backend

OperatorType = Literal[
    "ZERO",
    "ENERGY",
    "ELECTRIC_DIPOLE",
    "ELECTRIC_QUADRUPOLE",
    "ELECTRIC_OCTUPOLE",
    "MAGNETIC_DIPOLE",
    "DIAMAGNETIC",
    "ARBITRARY",
]

CPPOperatorType = _backend.OperatorType


def get_cpp_operator_type(operator_type: OperatorType) -> CPPOperatorType:
    if operator_type.upper() not in get_args(OperatorType):
        raise ValueError(f"Unknown operator_type '{operator_type}', should be one of {OperatorType}")
    try:
        return getattr(CPPOperatorType, f"{operator_type.upper()}")
    except AttributeError as err:
        raise ValueError(f"Unknown operator_type 'OperatorType.{operator_type.upper()}'") from err
