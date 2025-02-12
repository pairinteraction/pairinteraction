from collections.abc import Callable
from typing import TYPE_CHECKING, Any, Literal, get_args

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer, UnionCPPDiagonalizer
from pairinteraction.backend._wrapped.OperatorType import CPPOperatorType, OperatorType
from pairinteraction.backend._wrapped.Parity import CPPParity, Parity

if TYPE_CHECKING:
    from pairinteraction.backend._wrapped.system.System import System


def get_cpp_diagonalize(system: "System") -> Callable:
    try:
        return getattr(_backend, f"diagonalize{type(system).__name__}")
    except AttributeError as err:
        raise ValueError(f"Unknown diagonalize function for {type(system).__name__}") from err


def get_cpp_diagonalizer(diagonalizer: Diagonalizer, cpp_system: Any) -> UnionCPPDiagonalizer:
    if diagonalizer.capitalize() not in get_args(Diagonalizer):
        raise ValueError(f"Unknown diagonalizer '{diagonalizer}', should be one of {Diagonalizer}")
    type_ = get_type_of_system(cpp_system)
    try:
        return getattr(_backend, f"Diagonalizer{diagonalizer.capitalize()}{type_}")()
    except AttributeError as err:
        raise ValueError(f"Unknown diagonalizer 'Diagonalizer{diagonalizer.capitalize()}{type_}'") from err


def get_type_of_system(cpp_system: Any) -> Literal["Complex", "Real"]:
    if type(cpp_system).__name__.endswith("Complex"):
        return "Complex"
    if type(cpp_system).__name__.endswith("Real"):
        return "Real"
    raise ValueError("Unknown type of system")


def get_cpp_operator_type(operator_type: OperatorType) -> CPPOperatorType:
    if operator_type.upper() not in get_args(OperatorType):
        raise ValueError(f"Unknown operator_type '{operator_type}', should be one of {OperatorType}")
    try:
        return getattr(CPPOperatorType, f"{operator_type.upper()}")
    except AttributeError as err:
        raise ValueError(f"Unknown operator_type 'OperatorType.{operator_type.upper()}'") from err


def get_cpp_parity(parity: Parity) -> CPPParity:
    if parity.upper() not in get_args(Parity):
        raise ValueError(f"Unknown parity '{parity}', should be one of {Parity}")
    try:
        return getattr(CPPParity, f"{parity.upper()}")
    except AttributeError as err:
        raise ValueError(f"Unknown parity 'Parity.{parity.upper()}'") from err
