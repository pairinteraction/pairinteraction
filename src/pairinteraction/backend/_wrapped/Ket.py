from abc import ABC
from typing import TYPE_CHECKING, Union, overload

import pairinteraction.backend._backend as _backend
from pairinteraction.backend._wrapped.Parity import Parity
from pairinteraction.unit_system import Qty

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

UnionCPPKet = Union[_backend.KetFloat, _backend.KetDouble]


class KetBase(ABC):
    _cpp: UnionCPPKet

    @classmethod
    def _from_cpp_object(cls, cpp_obj: UnionCPPKet):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def __str__(self) -> str:
        return self.label

    @property
    def label(self) -> str:
        return self._cpp.get_label()

    @property
    def m(self) -> float:
        return self._cpp.get_quantum_number_m()

    @property
    def f(self) -> float:
        return self._cpp.get_quantum_number_f()

    @property
    def parity(self) -> Parity:
        parity_cpp = self._cpp.get_parity()
        return parity_cpp.name

    @property
    def energy(self) -> "PlainQuantity[float]":
        return self.get_energy()

    @overload
    def get_energy(self) -> "PlainQuantity[float]": ...

    @overload
    def get_energy(self, unit: str) -> float: ...

    def get_energy(self, unit: str = "pint"):
        energy_au = self._cpp.get_energy()
        energy = Qty.from_base(energy_au, "energy")
        return energy.to_unit(unit)
