# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any

from pairinteraction.ket import KetPair, KetPairReal
from pairinteraction.state.state import StateBase

if TYPE_CHECKING:
    from pairinteraction.basis import BasisPair, BasisPairReal


class StatePair(StateBase["BasisPair", KetPair]):
    """Pair state of two atoms.

    Currently StatePair objects don't offer any additional functionality.

    """

    _basis: "BasisPair"
    _ket_class = KetPair

    def get_amplitude(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_amplitude not implemented yet")

    def get_overlap(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_overlap not implemented yet")

    def get_matrix_element(self, other: Any, *args: Any, **kwargs: Any) -> Any:
        raise NotImplementedError("StatePair.get_matrix_element not implemented yet")


class StatePairReal(StatePair):
    _basis: "BasisPairReal"
    _ket_class = KetPairReal
