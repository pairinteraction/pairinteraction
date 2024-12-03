from abc import ABC
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar, Union, get_args, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.Parity import Parity
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    SelfKet_t = TypeVar("SelfKet_t", bound="Ket")

UnionCPPKet = Union[_backend.KetFloat, _backend.KetDouble]
UnionTypeCPPKetCreator = Any  # is supposed to be type[Ket(Atom|ClassicalLight)Creator(Float|Double)]


class KetBase(ABC):
    _cpp: UnionCPPKet
    _cpp_creator: ClassVar[UnionTypeCPPKetCreator]

    @classmethod
    def _from_cpp_object(cls: "type[SelfKet_t]", cpp_obj: UnionCPPKet) -> "SelfKet_t":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def __repr__(self) -> str:
        return self.label

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
        parity = parity_cpp.name
        if parity in get_args(Parity):
            return parity  # type: ignore [return-value]
        raise ValueError(f"Unknown parity {parity}")

    @property
    def energy(self):
        return self.get_energy()

    @overload
    def get_energy(self) -> "PlainQuantity[float]": ...

    @overload
    def get_energy(self, unit: str) -> float: ...

    def get_energy(self, unit: str = "pint"):
        energy_au = self._cpp.get_energy()
        energy = QuantityScalar.from_base(energy_au, "ENERGY")
        return energy.to_unit(unit)


Ket = KetBase
