from typing import TYPE_CHECKING, Any, Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.Ket import KetBase

if TYPE_CHECKING:
    pass

UnionCPPKetCombined = Union[
    _backend.KetCombinedFloat,
    _backend.KetCombinedComplexFloat,
    _backend.KetCombinedDouble,
    _backend.KetCombinedComplexDouble,
]
UnionTypeCPPKetCombinedCreator = Any


class KetCombinedBase(KetBase):
    _cpp: UnionCPPKetCombined
    _cpp_creator: UnionTypeCPPKetCombinedCreator


class KetCombinedFloat(KetCombinedBase):
    _cpp: _backend.KetCombinedFloat
    _cpp_creator = None


class KetCombinedComplexFloat(KetCombinedBase):
    _cpp: _backend.KetCombinedComplexFloat
    _cpp_creator = None


class KetCombinedDouble(KetCombinedBase):
    _cpp: _backend.KetCombinedDouble
    _cpp_creator = None


class KetCombinedComplexDouble(KetCombinedBase):
    _cpp: _backend.KetCombinedComplexDouble
    _cpp_creator = None
