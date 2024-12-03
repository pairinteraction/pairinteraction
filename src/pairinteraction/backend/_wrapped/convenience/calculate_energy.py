from typing import TYPE_CHECKING, Optional, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.ket.KetAtom import KetAtom
from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtom
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity


@overload
def calculate_energy(ket: KetAtom) -> "PlainQuantity[float]": ...


@overload
def calculate_energy(ket: KetAtom, system: SystemAtom) -> "PlainQuantity[float]": ...


@overload
def calculate_energy(ket: KetAtom, *, unit: str) -> float: ...


@overload
def calculate_energy(ket: KetAtom, system: SystemAtom, *, unit: str) -> float: ...


def calculate_energy(ket: KetAtom, system: Optional[SystemAtom] = None, *, unit: str = "pint"):
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
    return QuantityScalar.from_base(energy_au, "ENERGY").to_unit(unit)
