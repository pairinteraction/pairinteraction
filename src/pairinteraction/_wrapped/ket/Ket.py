from abc import ABC
from typing import TYPE_CHECKING, ClassVar, Literal, Optional, Union, get_args, overload

from pairinteraction import _backend
from pairinteraction._wrapped.cpp_types import Parity
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.units import PintFloat

UnionCPPKet = Union[_backend.KetAtom, _backend.KetPairComplex, _backend.KetPairReal]
UnionTypeCPPKetCreator = type[_backend.KetAtomCreator]  # since there currently is no KetPairCreator


class KetBase(ABC):
    """Base class for all Ket objects.

    The ket objects are meant to represent mathematically the canonical basis states, with respect to which
    the coefficient matrix of the basis objects are defined.
    For single atoms we simply choose the atomic states defined by their quantum numbers,
    therefore all KetAtom objects are orthogonal to each other.
    For pair systems, we choose the product states of the single atom eigenstates, which depends on the system
    and the applied fields. Thus for different pair systems the KetPair objects are not necessarily orthogonal anymore.

    All ket objects share a few common attributes and methods, that are defined in this base class.
    E.g. each ket has a total momentum quantum number f, a magnetic quantum number m, a parity, an energy,
    as well as a label that represents the ket.
    """

    _cpp: UnionCPPKet
    _cpp_creator: ClassVar[UnionTypeCPPKetCreator]

    @classmethod
    def _from_cpp_object(cls: "type[Self]", cpp_obj: UnionCPPKet) -> "Self":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def __repr__(self) -> str:
        return self.get_label("raw")

    def __str__(self) -> str:
        return self.get_label("ket")

    def get_label(self, format: Literal["ket", "bra", "raw"]) -> str:
        """Label representing the ket."""
        raw = self._cpp.get_label()
        if format == "raw":
            return raw
        if format == "ket":
            return f"|{raw}⟩"
        if format == "bra":
            return f"⟨{raw}|"
        raise ValueError(f"Unknown format {format}")

    @property
    def m(self) -> float:
        """The magnetic quantum number m (int or half-int)."""
        return self._cpp.get_quantum_number_m()

    @property
    def f(self) -> float:
        """The total momentum quantum number f (int or half-int)."""
        return self._cpp.get_quantum_number_f()

    @property
    def parity(self) -> Parity:
        """The parity of the ket."""
        parity_cpp = self._cpp.get_parity()
        parity = parity_cpp.name
        if parity in get_args(Parity):
            return parity  # type: ignore [return-value]
        raise ValueError(f"Unknown parity {parity}")

    @property
    def energy(self) -> "PintFloat":
        """The energy of the ket: E=I-Ry/nu^2."""
        return self.get_energy()

    @overload
    def get_energy(self, unit: None = None) -> "PintFloat": ...

    @overload
    def get_energy(self, unit: str) -> float: ...

    def get_energy(self, unit: Optional[str] = None) -> Union[float, "PintFloat"]:
        """Get the energy of the ket in the given unit.

        Args:
            unit: The unit to which to convert the energy to.
                Default None will return a `pint.Quantity`.

        Returns:
            The energy as float if a unit was given, otherwise a `pint.Quantity`.

        """
        energy_au = self._cpp.get_energy()
        energy = QuantityScalar.from_base_unit(energy_au, "ENERGY")
        return energy.to_pint_or_unit(unit)
