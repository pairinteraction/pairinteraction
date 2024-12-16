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
