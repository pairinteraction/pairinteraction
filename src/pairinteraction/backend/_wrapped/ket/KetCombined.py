from typing import Any, ClassVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.ket.Ket import KetBase

UnionCPPKetCombined = Union[
    _backend.KetCombinedFloat,
    _backend.KetCombinedComplexFloat,
    _backend.KetCombinedDouble,
    _backend.KetCombinedComplexDouble,
]
UnionTypeCPPKetCombinedCreator = Any


class KetCombinedBase(KetBase):
    _cpp: UnionCPPKetCombined  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator: ClassVar[UnionTypeCPPKetCombinedCreator]


class KetCombinedFloat(KetCombinedBase):
    _cpp: _backend.KetCombinedFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None


class KetCombinedComplexFloat(KetCombinedBase):
    _cpp: _backend.KetCombinedComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None


class KetCombinedDouble(KetCombinedBase):
    _cpp: _backend.KetCombinedDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None


class KetCombinedComplexDouble(KetCombinedBase):
    _cpp: _backend.KetCombinedComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = None
