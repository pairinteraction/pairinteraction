from collections.abc import Sequence
from typing import TYPE_CHECKING, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer
from pairinteraction.backend._wrapped.get_functions import get_cpp_diagonalize, get_cpp_diagonalizer, get_cpp_range
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    from pairinteraction.backend._wrapped.system.System import System


def diagonalize(
    systems: Sequence["System"],
    diagonalizer: Diagonalizer = "Eigen",
    precision: int = 12,
    sort_by_energy: bool = True,
    energy_range: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
    energy_unit: str = "pint",
) -> None:
    cpp_systems = [s._cpp for s in systems]  # type: ignore [reportPrivateUsage]
    cpp_diagonalize_fct = get_cpp_diagonalize(systems[0])

    cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, cpp_systems[0])
    if energy_range is None:
        cpp_diagonalize_fct(cpp_systems, cpp_diagonalizer, precision)
    else:
        min_energy_au = QuantityScalar(energy_range[0], energy_unit).to_base("ENERGY")
        max_energy_au = QuantityScalar(energy_range[1], energy_unit).to_base("ENERGY")
        cpp_range_class = get_cpp_range(cpp_systems[0])
        cpp_range = cpp_range_class(min_energy_au, max_energy_au)
        cpp_diagonalize_fct(cpp_systems, cpp_diagonalizer, precision, cpp_range)
    for system, cpp_system in zip(systems, cpp_systems):
        if sort_by_energy:
            sorter = cpp_system.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            cpp_system.transform(sorter)
        system._cpp = cpp_system  # type: ignore [reportPrivateUsage]
        system.update_basis()
