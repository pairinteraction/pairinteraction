from abc import ABC
from functools import cached_property
from typing import TYPE_CHECKING, Any, Generic, TypeVar

import scipy.sparse

from pairinteraction.backend._wrapped.Ket import KetBase

if TYPE_CHECKING:
    pass

# UnionCPPBasis = Union[
#     _backend.BasisBasis...,
# ]
UnionCPPBasis = Any

# UnionTypeCPPBasisCreator = Union[
#     type[_backend.BasisAtomCreatorFloat],
# ]
UnionTypeCPPBasisCreator = Any

# UnionTypeKet = Union[type[KetAtomFloat], type[KetAtomDouble]]
UnionTypeKet = Any

Ket_t = TypeVar("Ket_t", bound="KetBase")


class BasisBase(ABC, Generic[Ket_t]):
    _cpp: UnionCPPBasis
    _cpp_creator: UnionTypeCPPBasisCreator
    _TypeKet: UnionTypeKet

    @classmethod
    def _from_cpp_object(
        cls,
        cpp_obj: UnionCPPBasis,
    ):
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    @cached_property
    def kets(self) -> list[Ket_t]:
        """Return a list containing the kets of the basis."""
        kets = [self._TypeKet._from_cpp_object(ket) for ket in self._cpp.get_kets()]
        return kets

    @property
    def number_of_states(self) -> int:
        return self._cpp.get_number_of_states()

    @property
    def number_of_kets(self) -> int:
        return self._cpp.get_number_of_kets()

    @property
    def coefficients(self) -> scipy.sparse.csr_matrix:
        return self._cpp.get_coefficients()
