from typing import Optional

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.KetAtom import KetAtomBase
from pairinteraction.backend._wrapped.SystemAtom import SystemAtomBase
from pairinteraction.unit_system import Qty


def calculate_energy(ket: KetAtomBase, system: Optional[SystemAtomBase] = None, unit: str = "pint") -> float:
    """Calculate the energy of a ket state.

    Args:
        ket: The ket state.
        system (optional): The system to calculate the energy in. If no system is provided, simply return ket.get_energy

    Returns:
        The energy of the ket state.
    """
    if system is None:
        energy_au = _backend.calculate_energy(ket._cpp)
    else:
        energy_au = _backend.calculate_energy(ket._cpp, system._cpp)
    return Qty.from_base(energy_au, "energy").to_unit(unit)
