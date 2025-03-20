from abc import ABC
from functools import cached_property
from typing import TYPE_CHECKING, ClassVar, Generic, TypeVar, Union

import numpy as np

from pairinteraction import _backend
from pairinteraction._wrapped.ket.Ket import KetBase

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

KetType = TypeVar("KetType", bound=KetBase, covariant=True)
UnionCPPBasis = Union[
    _backend.BasisAtomReal, _backend.BasisAtomComplex, _backend.BasisPairReal, _backend.BasisPairComplex
]
UnionTypeCPPBasisCreator = Union[
    type[_backend.BasisAtomCreatorReal],
    type[_backend.BasisAtomCreatorComplex],
    type[_backend.BasisPairCreatorReal],
    type[_backend.BasisPairCreatorComplex],
]


class BasisBase(ABC, Generic[KetType]):
    """Base class for all Basis objects.

    The basis objects are meant to represent a set of kets, that span a Hilbert space and store a coefficient matrix,
    that describe the basis states in terms of the kets.

    All basis objects share a few common attributes and methods, that are defined in this base class, e.g.:
        - the number of kets and states,
        - the kets of the basis,
        - the coefficients stored as scipy sparse matrix,
        - ...
    """

    _cpp: UnionCPPBasis
    _cpp_creator: ClassVar[UnionTypeCPPBasisCreator]
    _TypeKet: type[KetType]  # should by ClassVar, but cannot be nested yet

    @classmethod
    def _from_cpp_object(cls: "type[Self]", cpp_obj: UnionCPPBasis) -> "Self":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def __repr__(self) -> str:
        args = f"{self.kets[0]} ... {self.kets[-1]}"
        return f"{type(self).__name__}({args})"

    def __str__(self) -> str:
        return self.__repr__().replace("Real", "").replace("Complex", "")

    @cached_property
    def kets(self) -> list[KetType]:
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
    def coefficients(self) -> "csr_matrix":
        return self._cpp.get_coefficients()

    def get_corresponding_state(self: "Self", ket_or_index: Union[KetType, int]) -> "Self":
        if isinstance(ket_or_index, (int, np.integer)):
            cpp_basis = self._cpp.get_corresponding_state(ket_or_index)
        else:
            cpp_basis = self._cpp.get_corresponding_state(ket_or_index._cpp)  # type: ignore [arg-type]
        return type(self)._from_cpp_object(cpp_basis)

    def get_corresponding_state_index(self, ket_or_index: Union[KetType, int]) -> int:
        if isinstance(ket_or_index, (int, np.integer)):
            return self._cpp.get_corresponding_state_index(ket_or_index)
        return self._cpp.get_corresponding_state_index(ket_or_index._cpp)  # type: ignore [arg-type]

    def get_corresponding_ket(self: "Self", state_or_index: Union["Self", int]) -> KetType:
        if isinstance(state_or_index, (int, np.integer)):
            cpp_ket = self._cpp.get_corresponding_ket(state_or_index)
        else:
            cpp_ket = self._cpp.get_corresponding_ket(state_or_index._cpp)  # type: ignore [arg-type]
        return self._TypeKet._from_cpp_object(cpp_ket)

    def get_corresponding_ket_index(self, state_or_index: Union["Self", int]) -> int:
        if isinstance(state_or_index, (int, np.integer)):
            return self._cpp.get_corresponding_ket_index(state_or_index)
        return self._cpp.get_corresponding_ket_index(state_or_index._cpp)  # type: ignore [arg-type]
