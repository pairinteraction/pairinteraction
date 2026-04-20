# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, TypeAlias, TypeGuard

import numpy as np

from pairinteraction import _backend
from pairinteraction.ket import KetPair, KetPairReal
from pairinteraction.state.state_atom import StateAtom
from pairinteraction.state.state_base import StateBase

if TYPE_CHECKING:
    from collections.abc import Sequence

    StatePairLike: TypeAlias = "StatePair | tuple[StateAtom, StateAtom] | Sequence[StateAtom]"

logger = logging.getLogger(__name__)


def is_state_pair_like(obj: Any) -> TypeGuard[StatePairLike]:
    return isinstance(obj, StatePair) or is_state_atom_tuple(obj)


def is_state_atom_tuple(obj: Any) -> TypeGuard[tuple[StateAtom, StateAtom]]:
    return hasattr(obj, "__len__") and len(obj) == 2 and all(isinstance(x, StateAtom) for x in obj)


class StatePair(StateBase[KetPair]):
    """Pair state of two atoms.

    Currently StatePair objects don't offer any additional functionality.

    """

    _cpp: _backend.BasisPairComplex
    _ket_class = KetPair

    def __init__(self, ket: KetPair, basis: Any) -> None:
        """Initialize a state object representing a ket in a given basis.

        Args:
            ket: The ket to represent in the state.
            basis: The basis to which the state belongs.

        """
        raise NotImplementedError(
            "StatePair objects cannot be created directly. "
            "You can use `basis_pair.get_corresponding_state(ket)` or `basis_pair.get_state(i)` instead."
        )

    def get_label(self, max_kets: int = 3) -> str:
        """Label representing the state.

        Args:
            max_kets: Maximum number of kets to include in the label.

        Returns:
            The label of the ket in the given format.

        """
        if isinstance(self._cpp, _backend.BasisPairComplex):
            from pairinteraction.basis import BasisAtom, BasisPair
            from pairinteraction.system import SystemAtom

            basis_atom_class = BasisAtom
            basis_pair_class = BasisPair
            basis_system_class = SystemAtom
        else:
            from pairinteraction.basis import BasisAtomReal, BasisPairReal
            from pairinteraction.system import SystemAtomReal

            basis_atom_class = BasisAtomReal
            basis_pair_class = BasisPairReal
            basis_system_class = SystemAtomReal

        basis_atoms_cpp = [self._cpp.get_basis1(), self._cpp.get_basis2()]
        basis_atoms = [
            basis_atom_class._from_cpp_object(basis_atom_cpp.get_canonical_basis())
            for basis_atom_cpp in basis_atoms_cpp
        ]
        system_atoms = [basis_system_class(basis_atom).diagonalize() for basis_atom in basis_atoms]

        # TODO this often does not work in the gui
        ket_pair = self.get_corresponding_ket()
        # TODO how to automatically choose the energy range?
        energy_range = (ket_pair.get_energy("GHz") - 10, ket_pair.get_energy("GHz") + 10)
        canonical_basis_pair = basis_pair_class(system_atoms, energy=energy_range, energy_unit="GHz")

        amplitudes = canonical_basis_pair.get_amplitudes(self)
        overlaps = canonical_basis_pair.get_overlaps(self)

        largest_inds = np.argpartition(overlaps, -max_kets)[-max_kets:]
        largest_inds = largest_inds[np.argsort(overlaps[largest_inds])[::-1]]

        if self.norm**2 - np.sum(overlaps) > overlaps[largest_inds[-1]]:
            logger.warning("The label may be inaccurate due to a to small canonical pair basis.")

        label = ""
        accumulated_ov = 0.0
        for ind in largest_inds:
            coeff = np.real_if_close(amplitudes[ind])
            canonical_ket = canonical_basis_pair.get_ket(ind)
            label += f"{coeff:.2f} |{canonical_ket.get_label()}⟩"
            accumulated_ov += overlaps[ind]
            label += " + "
            if accumulated_ov > (0.95 * self.norm**2):
                break

        if accumulated_ov <= self.norm**2 - 100 * np.finfo(float).eps:
            label += "..."
        else:
            label = label[:-3]  # Remove the last " + "

        return label.replace("+ -", "- ")

    def get_amplitude(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_amplitude not implemented yet")

    def get_overlap(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_overlap not implemented yet")

    def get_matrix_element(self, other: Any, *args: Any, **kwargs: Any) -> Any:
        raise NotImplementedError("StatePair.get_matrix_element not implemented yet")


class StatePairReal(StatePair):
    _cpp: _backend.BasisPairReal  # type: ignore [assignment]
    _ket_class = KetPairReal
