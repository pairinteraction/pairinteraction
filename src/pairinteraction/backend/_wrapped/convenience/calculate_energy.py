from typing import TYPE_CHECKING, Any, Optional, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.ket.KetAtom import KetAtomBase
from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtomBase
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity


@overload
def calculate_energy(ket: KetAtomBase, system: Optional[SystemAtomBase[Any]] = None) -> "PlainQuantity[float]": ...


@overload
def calculate_energy(ket: KetAtomBase, system: Optional[SystemAtomBase[Any]], unit: str) -> float: ...


def calculate_energy(ket: KetAtomBase, system: Optional[SystemAtomBase[Any]] = None, unit: str = "pint"):
    """Calculate the energy of a ket state.

    Args:
        ket: The ket state.
        system (optional): The system to calculate the energy in. If no system is provided, simply return ket.get_energy

    Returns:
        The energy of the ket state.
    """
    if system is None:
        energy_au = _backend.calculate_energy(ket._cpp)  # type: ignore [reportPrivateUsage]
    else:
        energy_au = _backend.calculate_energy(ket._cpp, system._cpp)  # type: ignore
    return Qty.from_base(energy_au, "energy").to_unit(unit)
