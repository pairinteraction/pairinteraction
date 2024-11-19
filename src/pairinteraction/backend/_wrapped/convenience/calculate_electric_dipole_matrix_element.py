from typing import Optional, TypeVar

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.ket.KetAtom import KetAtomDouble, KetAtomFloat
from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtomBase

Ket_t = TypeVar("Ket_t", KetAtomDouble, KetAtomFloat)


def calculate_electric_dipole_matrix_element(
    initial_ket: Ket_t, final_ket: Ket_t, q: int, system: Optional[SystemAtomBase] = None
) -> float:
    """Calculate the electric dipole matrix element between two kets.

    Args:
        initial_ket: The initial ket state.
        final_ket: The final ket state.
        q: The dipole matrix element of the d^q operator.
        system (optional): The system to calculate the matrix element in.
            If no system is provided, simply return the matrix element between the two kets.

    Returns:
        The electric dipole matrix element between the two kets.
    """
    if system is None:
        return _backend.calculate_electric_dipole_matrix_element(initial_ket._cpp, final_ket._cpp, q)  # type: ignore

    return _backend.calculate_electric_dipole_matrix_element(initial_ket._cpp, final_ket._cpp, system._cpp, q)  # type: ignore
