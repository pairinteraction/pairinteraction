from typing import Any, Literal, Union, get_args

import pairinteraction.backend._backend as _backend

Diagonalizer = Literal["Eigen", "Lapacke", "Feast"]

UnionCPPDiagonalizer = Union[
    _backend.DiagonalizerInterfaceFloat,
    _backend.DiagonalizerInterfaceComplexFloat,
    _backend.DiagonalizerInterfaceDouble,
    _backend.DiagonalizerInterfaceComplexDouble,
]


def get_cpp_diagonalizer(diagonalizer: Diagonalizer, cpp_system: Any) -> UnionCPPDiagonalizer:
    if isinstance(diagonalizer, UnionCPPDiagonalizer):
        return diagonalizer
    if diagonalizer.capitalize() not in get_args(Diagonalizer):
        raise ValueError(f"Unknown diagonalizer '{diagonalizer}', should be one of {Diagonalizer}")
    type_ = get_type_of_system(cpp_system)
    try:
        return getattr(_backend, f"Diagonalizer{diagonalizer.capitalize()}{type_}")()
    except AttributeError as err:
        raise ValueError(f"Unknown diagonalizer 'Diagonalizer{diagonalizer.capitalize()}{type_}'") from err


def get_type_of_system(cpp_system: Any) -> Literal["ComplexFloat", "Float", "ComplexDouble", "Double"]:
    if type(cpp_system).__name__.endswith("ComplexFloat"):
        return "ComplexFloat"
    if type(cpp_system).__name__.endswith("Float"):
        return "Float"
    if type(cpp_system).__name__.endswith("ComplexDouble"):
        return "ComplexDouble"
    if type(cpp_system).__name__.endswith("Double"):
        return "Double"
    raise ValueError("Unknown type of system")
