from collections.abc import Callable
from typing import TYPE_CHECKING, Any, Literal, Optional, Union, get_args

from pairinteraction import _backend

if TYPE_CHECKING:
    from pairinteraction._wrapped.system.System import System


Diagonalizer = Literal["eigen", "lapacke_evd", "lapacke_evr", "feast"]
FloatType = Literal["float32", "float64"]
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
Parity = Literal["EVEN", "ODD", "UNKNOWN"]


UnionCPPDiagonalizer = Union[_backend.DiagonalizerInterfaceReal, _backend.DiagonalizerInterfaceComplex]
CPPFloatType = _backend.FloatType
CPPOperatorType = _backend.OperatorType
CPPParity = _backend.Parity


def get_cpp_diagonalize(system: "System") -> Callable:
    try:
        return getattr(_backend, f"diagonalize{type(system).__name__}")
    except AttributeError as err:
        raise ValueError(f"Unknown diagonalize function for {type(system).__name__}") from err


def get_cpp_diagonalizer(
    diagonalizer: Diagonalizer,
    cpp_system: Any,
    float_type: FloatType,
    m0: Optional[int] = None,
) -> UnionCPPDiagonalizer:
    if diagonalizer not in get_args(Diagonalizer):
        raise ValueError(f"Unknown diagonalizer '{diagonalizer}', should be one of {Diagonalizer}")
    if diagonalizer == "feast" and m0 is None:
        raise ValueError("m0 must be specified for the 'feast' diagonalizer")
    elif diagonalizer != "feast" and m0 is not None:
        raise ValueError("m0 must not be specified if the diagonalizer is not 'feast'")

    type_ = get_type_of_system(cpp_system)
    diagonalizer_ = "".join([s.capitalize() for s in diagonalizer.split("_")])
    try:
        diagonalizer_class = getattr(_backend, f"Diagonalizer{diagonalizer_}{type_}")
    except AttributeError as err:
        raise ValueError(f"Unknown diagonalizer 'Diagonalizer{diagonalizer_}{type_}'") from err

    cpp_float_type = get_cpp_float_type(float_type)
    if diagonalizer == "feast":
        return diagonalizer_class(m0=m0, float_type=cpp_float_type)
    return diagonalizer_class(float_type=cpp_float_type)


def get_type_of_system(cpp_system: Any) -> Literal["Complex", "Real"]:
    if type(cpp_system).__name__.endswith("Complex"):
        return "Complex"
    if type(cpp_system).__name__.endswith("Real"):
        return "Real"
    raise ValueError("Unknown type of system")


def get_cpp_operator_type(operator_type: OperatorType) -> CPPOperatorType:
    operator_type = operator_type.upper()
    if operator_type not in get_args(OperatorType):
        raise ValueError(f"Unknown operator_type '{operator_type}', should be one of {OperatorType}")
    try:
        return getattr(CPPOperatorType, f"{operator_type.upper()}")
    except AttributeError as err:
        raise ValueError(f"Unknown operator_type 'OperatorType.{operator_type.upper()}'") from err


def get_cpp_parity(parity: Parity) -> CPPParity:
    parity = parity.upper()
    if parity not in get_args(Parity):
        raise ValueError(f"Unknown parity '{parity}', should be one of {Parity}")
    try:
        return getattr(CPPParity, f"{parity.upper()}")
    except AttributeError as err:
        raise ValueError(f"Unknown parity 'Parity.{parity.upper()}'") from err


def get_cpp_float_type(float_type: FloatType) -> CPPFloatType:
    float_type = float_type.lower()
    if float_type not in get_args(FloatType):
        raise ValueError(f"Unknown float_type '{float_type}', should be one of {FloatType}")
    try:
        return getattr(CPPFloatType, f"{float_type.upper()}")
    except AttributeError as err:
        raise ValueError(f"Unknown float_type 'FloatType.{float_type.upper()}'") from err
