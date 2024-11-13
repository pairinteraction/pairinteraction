from abc import ABC
from typing import Union

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.BasisAtom import BasisAtomBase


class SystemAtomBase(ABC):
    _CPPSystemAtom: Union[
        type[_backend.SystemAtomFloat],
        type[_backend.SystemAtomComplexFloat],
        type[_backend.SystemAtomDouble],
        type[_backend.SystemAtomComplexDouble],
    ]

    def __init__(self, basis: "BasisAtomBase") -> None:
        self.basis = basis
        self._cpp = self._CPPSystemAtom(basis._cpp)


class SystemAtomFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomFloat
    _CPPSystemAtom = _backend.SystemAtomFloat


class SystemAtomComplexFloat(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexFloat
    _CPPSystemAtom = _backend.SystemAtomComplexFloat


class SystemAtomDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomDouble
    _CPPSystemAtom = _backend.SystemAtomDouble


class SystemAtomComplexDouble(SystemAtomBase):
    _cpp: _backend.SystemAtomComplexDouble
    _CPPSystemAtom = _backend.SystemAtomComplexDouble
