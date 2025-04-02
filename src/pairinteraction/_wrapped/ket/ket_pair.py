# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Sequence
from typing import Any, ClassVar, Union

from typing_extensions import TypeGuard

from pairinteraction import _backend
from pairinteraction._wrapped.ket.ket import KetBase
from pairinteraction._wrapped.ket.ket_atom import KetAtom

UnionCPPKetPair = Union[_backend.KetPairReal, _backend.KetPairComplex]
UnionTypeCPPKetPairCreator = Any

KetPairLike = Union["KetPairReal", "KetPairComplex", tuple["KetAtom", "KetAtom"], Sequence["KetAtom"]]


def is_ket_pair_like(obj: Any) -> TypeGuard[KetPairLike]:
    return isinstance(obj, KetPair) or (len(obj) == 2 and all(isinstance(x, KetAtom) for x in obj))


def is_ket_atom_tuple(obj: Any) -> TypeGuard[tuple["KetAtom", "KetAtom"]]:
    return len(obj) == 2 and all(isinstance(x, KetAtom) for x in obj)


class KetPair(KetBase):
    """Ket for a pair state of two atoms.

    For pair systems, we choose KetPair object as the product states of the single atom eigenstates.
    Thus, the Ket pair objects depend on the system and the applied fields.
    Therefore for different pair systems the KetPair objects are not necessarily orthogonal anymore.

    Currently one cannot create a KetPair object directly, but they are used in the background when creating a
    :class:`pairinteraction._wrapped.BasisPair` object.

    """

    _cpp: UnionCPPKetPair
    _cpp_creator: ClassVar[UnionTypeCPPKetPairCreator] = None

    def __init__(self) -> None:
        """Creating a KetPair object directly is not possible."""  # noqa: D401
        raise NotImplementedError("KetPair objects cannot be created directly.")


class KetPairReal(KetPair):
    _cpp: _backend.KetPairReal


class KetPairComplex(KetPair):
    _cpp: _backend.KetPairComplex
