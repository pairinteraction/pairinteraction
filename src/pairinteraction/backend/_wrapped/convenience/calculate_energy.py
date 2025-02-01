import warnings
from typing import TYPE_CHECKING, Optional, overload

import numpy as np

from pairinteraction.backend._wrapped.ket.KetAtom import KetAtom
from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtom

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
        unit (optional): The unit, in which to return the energy in. If not provided (or "pint"), return a pint.Quantity

    Returns:
        The energy of the ket state.

    """
    if system is None:
        return ket.get_energy(unit=unit)

    overlaps = system.get_eigenbasis().get_overlaps(ket)
    idx = np.argmax(overlaps)
    if overlaps[idx] < 0.5:
        warnings.warn(
            "The provided ket states does not correspond to an eigenstate of the system in a unique way.", stacklevel=2
        )

    return system.get_eigenvalues(unit=unit)[idx]
