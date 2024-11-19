from collections.abc import Sequence
from typing import TypeVar, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer, get_cpp_diagonalizer
from pairinteraction.backend._wrapped.system.SystemAtom import (
    SystemAtomBase,
)
from pairinteraction.backend._wrapped.system.SystemCombined import (
    SystemCombinedBase,
)

System_t = TypeVar("System_t", SystemAtomBase, SystemCombinedBase)


def diagonalize(
    systems: Sequence[System_t],
    diagonalizer: Diagonalizer = "Eigen",
    precision: int = 12,
    eigenvalue_range: Union[_backend.RangeFloat, _backend.RangeDouble, None] = None,
) -> list[System_t]:
    cpp_systems = [s._cpp for s in systems]
    diagonalize_fct = get_diagonalize_fct(systems[0])

    cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, cpp_systems[0])
    if eigenvalue_range is None:
        diagonalize_fct(cpp_systems, cpp_diagonalizer, precision)
    else:
        diagonalize_fct(cpp_systems, cpp_diagonalizer, precision, eigenvalue_range)
    for i, system in enumerate(systems):
        system._cpp = cpp_systems[i]
        system.update_basis()
    return list(systems)


def get_diagonalize_fct(system):
    try:
        return getattr(_backend, f"diagonalize{type(system).__name__}")
    except AttributeError as err:
        raise ValueError(f"Unknown diagonalize function for {type(system).__name__}") from err
