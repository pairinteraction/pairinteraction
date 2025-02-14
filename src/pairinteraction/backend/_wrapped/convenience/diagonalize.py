from collections.abc import Sequence
from typing import TYPE_CHECKING, Optional, TypeVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.Diagonalizer import Diagonalizer
from pairinteraction.backend._wrapped.get_functions import get_cpp_diagonalize, get_cpp_diagonalizer
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    from pairinteraction.backend._wrapped.system.System import System

    Quantity = TypeVar("Quantity", float, PlainQuantity[float])


def diagonalize(
    systems: Sequence["System"],
    diagonalizer: Diagonalizer = "Eigen",
    atol: float = 1e-6,
    sort_by_energy: bool = True,
    energy_range: tuple[Union["Quantity", None], Union["Quantity", None]] = (None, None),
    energy_unit: Optional[str] = None,
) -> None:
    cpp_systems = [s._cpp for s in systems]  # type: ignore [reportPrivateUsage]
    cpp_diagonalize_fct = get_cpp_diagonalize(systems[0])

    cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, cpp_systems[0])
    min_energy_au, max_energy_au = energy_range
    if min_energy_au is not None:
        min_energy_au = QuantityScalar(min_energy_au, energy_unit).to_base("ENERGY")
    if max_energy_au is not None:
        max_energy_au = QuantityScalar(max_energy_au, energy_unit).to_base("ENERGY")
    cpp_diagonalize_fct(cpp_systems, cpp_diagonalizer, min_energy_au, max_energy_au, atol)
    for system, cpp_system in zip(systems, cpp_systems):
        if sort_by_energy:
            sorter = cpp_system.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            cpp_system.transform(sorter)
        system._cpp = cpp_system  # type: ignore [reportPrivateUsage]
        system.update_basis()
