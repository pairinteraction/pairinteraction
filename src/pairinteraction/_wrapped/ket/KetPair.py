from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, ClassVar, Union

from pairinteraction import _backend
from pairinteraction._wrapped.ket.Ket import KetBase

if TYPE_CHECKING:
    from pairinteraction._wrapped.ket.KetAtom import KetAtom

UnionCPPKetPair = Union[_backend.KetPairReal, _backend.KetPairComplex]
UnionTypeCPPKetPairCreator = Any

KetPairLike = Union["KetPair", tuple["KetAtom", "KetAtom"], Sequence["KetAtom"]]


class KetPair(KetBase):
    """Ket for a pair state of two atoms.

    For pair systems, we choose KetPair object as the product states of the single atom eigenstates.
    Thus, the Ket pair objects depend on the system and the applied fields.
    Therefore for different pair systems the KetPair objects are not necessarily orthogonal anymore.

    Currently one cannot create a KetPair object directly, but they are used in the background when creating a
    :class:`pairinteraction._wrapped.BasisPair` object.

    """

    _cpp: UnionCPPKetPair
    _cpp_creator: ClassVar[UnionTypeCPPKetPairCreator]

    def __init__(self) -> None:
        """Creating a KetPair object directly is not possible."""  # noqa: D401
        raise NotImplementedError("KetPair objects cannot be created directly.")


class KetPairReal(KetPair):
    _cpp: _backend.KetPairReal
    _cpp_creator = None


class KetPairComplex(KetPair):
    _cpp: _backend.KetPairComplex
    _cpp_creator = None
