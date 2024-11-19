from collections.abc import Sequence
from typing import Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer, get_cpp_diagonalizer
from pairinteraction.backend._wrapped.system.System import UnionSystem


def diagonalize(
    systems: Sequence[UnionSystem],
    diagonalizer: Diagonalizer = "Eigen",
    precision: int = 12,
    eigenvalue_range: Union[_backend.RangeFloat, _backend.RangeDouble, None] = None,
    sort_by_energy: bool = True,
) -> None:
    cpp_systems = [s._cpp for s in systems]  # type: ignore [reportPrivateUsage]
    diagonalize_fct = get_cpp_diagonalize(systems[0])

    cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, cpp_systems[0])
    if eigenvalue_range is None:
        diagonalize_fct(cpp_systems, cpp_diagonalizer, precision)
    else:
        diagonalize_fct(cpp_systems, cpp_diagonalizer, precision, eigenvalue_range)
    for system, cpp_system in zip(systems, cpp_systems):
        if sort_by_energy:
            sorter = cpp_system.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            cpp_system.transform(sorter)
        system._cpp = cpp_system  # type: ignore [reportPrivateUsage]
        system.update_basis()


def get_cpp_diagonalize(system: UnionSystem):
    try:
        return getattr(_backend, f"diagonalize{type(system).__name__}")
    except AttributeError as err:
        raise ValueError(f"Unknown diagonalize function for {type(system).__name__}") from err
