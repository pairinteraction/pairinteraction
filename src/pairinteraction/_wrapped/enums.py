from typing import Literal, get_args

from pairinteraction._backend import (
    FloatType as CPPFloatType,
    OperatorType as CPPOperatorType,
    Parity as CPPParity,
)

FloatType = Literal["float32", "float64"]
OperatorType = Literal[
    "zero",
    "energy",
    "electric_dipole",
    "electric_quadrupole",
    "electric_quadrupole_zero",
    "electric_octupole",
    "magnetic_dipole",
    "arbitrary",
]
Parity = Literal["even", "odd", "unknown"]


def get_cpp_float_type(float_type: FloatType) -> CPPFloatType:
    if float_type not in get_args(FloatType):
        raise ValueError(f"Unknown float_type '{float_type}', should be one of {FloatType}")
    try:
        return getattr(CPPFloatType, f"{float_type.upper()}")  # type: ignore [no-any-return]
    except AttributeError as err:
        raise ValueError(f"Unknown float_type 'FloatType.{float_type.upper()}'") from err


def get_cpp_operator_type(operator_type: OperatorType) -> CPPOperatorType:
    if operator_type not in get_args(OperatorType):
        raise ValueError(f"Unknown operator_type '{operator_type}', should be one of {OperatorType}")
    try:
        return getattr(CPPOperatorType, f"{operator_type.upper()}")  # type: ignore [no-any-return]
    except AttributeError as err:
        raise ValueError(f"Unknown operator_type 'OperatorType.{operator_type.upper()}'") from err


def get_cpp_parity(parity: Parity) -> CPPParity:
    if parity not in get_args(Parity):
        raise ValueError(f"Unknown parity '{parity}', should be one of {Parity}")
    try:
        return getattr(CPPParity, f"{parity.upper()}")  # type: ignore [no-any-return]
    except AttributeError as err:
        raise ValueError(f"Unknown parity 'Parity.{parity.upper()}'") from err
