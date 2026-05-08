# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, Any, TypeAlias, TypeGuard

from pairinteraction.ket import KetPair, KetPairReal
from pairinteraction.state.state_atom import StateAtom
from pairinteraction.state.state_base import StateBase

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pairinteraction import _backend

StatePairLike: TypeAlias = "StatePair | tuple[StateAtom, StateAtom] | Sequence[StateAtom]"


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

    def get_amplitude(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_amplitude not implemented yet")

    def get_overlap(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_overlap not implemented yet")

    def get_matrix_element(self, other: Any, *args: Any, **kwargs: Any) -> Any:
        raise NotImplementedError("StatePair.get_matrix_element not implemented yet")


class StatePairReal(StatePair):
    _cpp: _backend.BasisPairReal  # type: ignore [assignment]
    _ket_class = KetPairReal
