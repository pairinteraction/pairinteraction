# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, TypeVar

from pairinteraction._wrapped.ket.ket_pair import KetPair, KetPairComplex, KetPairReal
from pairinteraction._wrapped.state.state import StateBase

if TYPE_CHECKING:
    from pairinteraction._wrapped.basis.basis_pair import BasisPair, BasisPairComplex, BasisPairReal

KetPairType = TypeVar("KetPairType", bound="KetPair", covariant=True)
BasisPairType = TypeVar("BasisPairType", bound="BasisPair[Any, Any]", covariant=True)


class StatePair(StateBase[BasisPairType, KetPairType]):
    """Pair state of two atoms.

    Currently StatePair objects don't offer any additional functionality.

    """

    def get_amplitude(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_amplitude not implemented yet")

    def get_overlap(self, other: Any) -> Any:
        raise NotImplementedError("StatePair.get_overlap not implemented yet")

    def get_matrix_element(self, other: Any, *args: Any, **kwargs: Any) -> Any:
        raise NotImplementedError("StatePair.get_matrix_element not implemented yet")


class StatePairReal(StatePair["BasisPairReal", "KetPairReal"]):
    _basis: "BasisPairReal"
    _TypeKet = KetPairReal


class StatePairComplex(StatePair["BasisPairComplex", "KetPairComplex"]):
    _basis: "BasisPairComplex"
    _TypeKet = KetPairComplex
