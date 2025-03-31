from typing import Any, Literal, Optional, Union, get_args

from pairinteraction import _backend
from pairinteraction._wrapped.enums import FloatType, get_cpp_float_type

Diagonalizer = Literal["eigen", "lapacke_evd", "lapacke_evr", "feast"]
UnionCPPDiagonalizer = Union[_backend.DiagonalizerInterfaceReal, _backend.DiagonalizerInterfaceComplex]


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
        return diagonalizer_class(m0=m0, float_type=cpp_float_type)  # type: ignore [no-any-return]
    return diagonalizer_class(float_type=cpp_float_type)  # type: ignore [no-any-return]


def get_type_of_system(cpp_system: Any) -> Literal["Complex", "Real"]:
    if type(cpp_system).__name__.endswith("Complex"):
        return "Complex"
    if type(cpp_system).__name__.endswith("Real"):
        return "Real"
    raise ValueError("Unknown type of system")
