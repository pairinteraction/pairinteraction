from abc import ABC
from functools import cached_property
from typing import TYPE_CHECKING, Any, ClassVar, Generic, TypeVar, Union, overload

from pairinteraction.backend._wrapped.ket.Ket import KetBase

if TYPE_CHECKING:
    import scipy.sparse
    from numpy.typing import ArrayLike

    SelfBasis_t = TypeVar("SelfBasis_t", bound="BasisBase[Any]")

Ket_t = TypeVar("Ket_t", bound=KetBase)
UnionCPPBasis = Any
# Basis(|Basis)(Atom|ClassicalLight|Combined)(Float|Double|ComplexFloat|ComplexDouble)
UnionTypeCPPBasisCreator = Any
# type[Basis(Atom|ClassicalLight|Combined)Creator(Float|Double|ComplexFloat|ComplexDouble)]


class BasisBase(ABC, Generic[Ket_t]):
    _cpp: UnionCPPBasis
    _cpp_creator: ClassVar[UnionTypeCPPBasisCreator]
    _TypeKet: type[Ket_t]  # should by ClassVar, but cannot be nested yet

    @classmethod
    def _from_cpp_object(cls: "type[SelfBasis_t]", cpp_obj: UnionCPPBasis) -> "SelfBasis_t":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    @cached_property
    def kets(self) -> list[Ket_t]:
        """Return a list containing the kets of the basis."""
        kets = [self._TypeKet._from_cpp_object(ket) for ket in self._cpp.get_kets()]  # type: ignore [reportPrivateUsage]
        return kets

    @property
    def number_of_states(self) -> int:
        return self._cpp.get_number_of_states()

    @property
    def number_of_kets(self) -> int:
        return self._cpp.get_number_of_kets()

    @property
    def coefficients(self) -> "scipy.sparse.csr_matrix":
        return self._cpp.get_coefficients()

    @overload
    def get_overlaps(self, ket_or_basis: Ket_t) -> "ArrayLike": ...

    @overload
    def get_overlaps(self, ket_or_basis: "BasisBase[Ket_t]") -> "scipy.sparse.csr_matrix": ...

    def get_overlaps(self, ket_or_basis: Union[Ket_t, "BasisBase[Ket_t]"]):
        return self._cpp.get_overlaps(ket_or_basis._cpp)

    def get_corresponding_state(self: "SelfBasis_t", ket_or_index: Union[KetBase, int]) -> "SelfBasis_t":
        if isinstance(ket_or_index, int):
            cpp_basis = self._from_cpp_object(self._cpp.get_corresponding_state(ket_or_index))
        else:
            cpp_basis = self._cpp.get_corresponding_state(ket_or_index._cpp)  # type: ignore [reportPrivateUsage]
        return type(self)._from_cpp_object(cpp_basis)

    def get_corresponding_state_index(self, ket_or_index: Union[KetBase, int]) -> int:
        if isinstance(ket_or_index, int):
            return self._cpp.get_corresponding_state_index(ket_or_index)
        return self._cpp.get_corresponding_state_index(ket_or_index._cpp)  # type: ignore [reportPrivateUsage]
