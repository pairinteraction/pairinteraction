from typing import TYPE_CHECKING, Optional, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.KetAtom import KetAtomBase
from pairinteraction.backend._wrapped.SystemAtom import SystemAtomBase
from pairinteraction.unit_system import convert_base_to_quantity

if TYPE_CHECKING:
    from pint.facets.plain import PlainUnit


def calculate_energy(
    ket: KetAtomBase, system: Optional[SystemAtomBase] = None, unit: Union[str, "PlainUnit"] = "GHz"
) -> float:
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
    return convert_base_to_quantity(energy_au, "energy", unit).magnitude
