from typing import Any, ClassVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.ket.Ket import KetBase

UnionCPPKetPair = Union[
    _backend.KetPairFloat,
    _backend.KetPairComplexFloat,
    _backend.KetPairDouble,
    _backend.KetPairComplexDouble,
]
UnionTypeCPPKetPairCreator = Any


class KetPairBase(KetBase):
    _cpp: UnionCPPKetPair  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator: ClassVar[UnionTypeCPPKetPairCreator]

    def __init__(self) -> None:
        """KetPair class.

        For pair systems, we choose KetPair object as the product states of the single atom eigenstates.
        Thus, the Ket pair objects depend on the system and the applied fields.
        Therefore for different pair systems the KetPair objects are not necessarily orthogonal anymore.

        Currently one cannot create a KetPair object directly, but they are used in the background when creating a
        :class:`pairinteraction.backend.double.BasisPair` object.
        """
        raise NotImplementedError("KetPair objects cannot be created directly.")


class KetPairFloat(KetPairBase):
    _cpp: _backend.KetPairFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None


class KetPairComplexFloat(KetPairBase):
    _cpp: _backend.KetPairComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None


class KetPairDouble(KetPairBase):
    _cpp: _backend.KetPairDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None


class KetPairComplexDouble(KetPairBase):
    _cpp: _backend.KetPairComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None


KetPair = KetPairBase
