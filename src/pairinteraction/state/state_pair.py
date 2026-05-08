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

    def get_label(
        self,
        stop_after_num_kets: int = 3,
        stop_after_accumulated_overlap: float = 0.95,
        considered_num_kets: int | None = None,
    ) -> str:
        """Label representing the state.

        Args:
            stop_after_num_kets: Maximum number of kets to include in the label.
            stop_after_accumulated_overlap: Stop including kets in the label,
                if the accumulated overlap of the included kets exceeds this value.
            considered_num_kets: The number of kets to consider in the pair basis.
                Default None uses a heuristic to determine a suitable number of kets.

        Returns:
            The label of the ket in the given format.

        """
        if isinstance(self._cpp, _backend.BasisPairComplex):
            from pairinteraction.basis import BasisAtom, BasisPair
            from pairinteraction.system import SystemAtom

            basis_atom_class = BasisAtom
            system_atom_class = SystemAtom
            basis_pair_class = BasisPair
        else:
            from pairinteraction.basis import BasisAtomReal, BasisPairReal
            from pairinteraction.system import SystemAtomReal

            basis_atom_class = BasisAtomReal
            system_atom_class = SystemAtomReal
            basis_pair_class = BasisPairReal

        basis_atoms_cpp = [self._cpp.get_basis1(), self._cpp.get_basis2()]
        basis_atoms = [
            basis_atom_class._from_cpp_object(basis_atom_cpp.canonicalized()) for basis_atom_cpp in basis_atoms_cpp
        ]
        system_atoms = [system_atom_class(basis_atom) for basis_atom in basis_atoms]

        coeffs = np.abs(self.get_coefficients())
        ket_pair = self.get_ket(int(np.argmax(coeffs)))
        ket_atom_tuple = tuple(state_atom.get_corresponding_ket() for state_atom in ket_pair.state_atoms)

        # heuristic to quickly find a basis, which includes the stop_after_num_kets most contributing kets
        is_converged = False
        considered_num_kets_list = [100, 1_000, 10_000] if considered_num_kets is None else [considered_num_kets]
        for number_of_kets in considered_num_kets_list:
            canonical_basis_pair = basis_pair_class.from_kets(
                ket_atom_tuple, system_atoms, number_of_kets=number_of_kets, warn_number_of_kets=False
            )

            amplitudes = canonical_basis_pair.get_matrix_elements(self, ("identity", "identity"), (0, 0)).m
            overlaps = np.abs(amplitudes) ** 2

            _stop_after_num_kets = min(stop_after_num_kets, len(overlaps))
            largest_inds = np.argpartition(overlaps, -_stop_after_num_kets)[-_stop_after_num_kets:]
            largest_inds = largest_inds[np.argsort(overlaps[largest_inds])[::-1]]

            acc_overlaps = np.cumsum(overlaps[largest_inds])
            max_ind_to_include = np.searchsorted(acc_overlaps, stop_after_accumulated_overlap * self.norm**2)
            largest_inds = largest_inds[: max_ind_to_include + 1]

            # if the overlap of the smallest ket we still include is larger than the remaining contributions,
            # we can be sure, that the label is accurate and the heuristic was successful.
            remaining_contributions = self.norm**2 - np.sum(overlaps)
            if overlaps[largest_inds[-1]] >= remaining_contributions:
                is_converged = True
                break

        label = ""
        accumulated_ov = 0.0
        for ind in largest_inds:
            coeff = np.real_if_close(amplitudes[ind])
            canonical_ket = canonical_basis_pair.get_ket(ind)
            label += f"{coeff:.2f} |{canonical_ket.get_label()}⟩"
            accumulated_ov += overlaps[ind]
            label += " + "

        if accumulated_ov <= self.norm**2 - 100 * np.finfo(float).eps:
            label += "..."
        else:
            label = label[:-3]  # Remove the last " + "
        label = label.replace("+ -", "- ")

        if not is_converged:
            logger.warning(
                "The label '%s' may not be printing the largest contributions. "
                "Consider calling get_label with a larger 'considered_num_kets'.",
                label,
            )

        return label

    def get_amplitude(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_amplitude not implemented yet")

    def get_overlap(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_overlap not implemented yet")

    def get_matrix_element(self, other: Any, *args: Any, **kwargs: Any) -> Any:
        raise NotImplementedError("StatePair.get_matrix_element not implemented yet")


class StatePairReal(StatePair):
    _cpp: _backend.BasisPairReal  # type: ignore [assignment]
    _ket_class = KetPairReal
