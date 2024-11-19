from typing import Optional

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.ket.KetAtom import KetAtomBase
from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtomBase


def calculate_electric_dipole_matrix_element(
    initial_ket: KetAtomBase, final_ket: KetAtomBase, q: int, system: Optional[SystemAtomBase] = None
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
        return _backend.calculate_electric_dipole_matrix_element(initial_ket._cpp, final_ket._cpp, q)

    return _backend.calculate_electric_dipole_matrix_element(initial_ket._cpp, final_ket._cpp, system._cpp, q)
