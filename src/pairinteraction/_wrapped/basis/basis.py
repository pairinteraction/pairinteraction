# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, ClassVar, Generic, TypeVar, Union

from pairinteraction import _backend

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction._wrapped.ket.ket import KetBase
    from pairinteraction._wrapped.state.state import StateBase

KetType = TypeVar("KetType", bound="KetBase", covariant=True)
StateType = TypeVar("StateType", bound="StateBase[Any, Any]", covariant=True)
UnionCPPBasis = Union[
    _backend.BasisAtomReal, _backend.BasisAtomComplex, _backend.BasisPairReal, _backend.BasisPairComplex
]
UnionTypeCPPBasisCreator = Union[
    type[_backend.BasisAtomCreatorReal],
    type[_backend.BasisAtomCreatorComplex],
    type[_backend.BasisPairCreatorReal],
    type[_backend.BasisPairCreatorComplex],
]


class BasisBase(ABC, Generic[KetType, StateType]):
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
    _TypeKet: type[KetType]  # should be ClassVar, but cannot be nested yet
    _TypeState: type[StateType]  # should be ClassVar, but cannot be nested yet

    @abstractmethod
    def __init__(self, *args: Any, **kwargs: Any) -> None: ...

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

    @property
    def kets(self) -> list[KetType]:
        """Return a list containing the kets of the basis."""
        return [self._TypeKet._from_cpp_object(ket) for ket in self._cpp.get_kets()]

    @property
    def states(self) -> list[StateType]:
        """Return a list containing the states of the basis."""
        states: list[StateType] = []
        for i in range(self.number_of_states):
            state_cpp = self._cpp.get_state(i)
            state_basis = type(self)._from_cpp_object(state_cpp)
            states.append(self._TypeState._from_basis_object(state_basis))
        return states

    @property
    def number_of_kets(self) -> int:
        """Return the number of kets in the basis."""
        return self._cpp.get_number_of_kets()

    @property
    def number_of_states(self) -> int:
        """Return the number of states in the basis."""
        return self._cpp.get_number_of_states()

    def get_coefficients(self) -> "csr_matrix":
        """Return the coefficients of the basis as a sparse matrix.

        The coefficients are stored in a sparse matrix with shape (number_of_kets, number_of_states),
        where the first index correspond to the kets and the second index correspond to the states.
        For example `basis.get_coefficients()[i, j]` is the i-th coefficient
        (i.e. the coefficient corresponding to the i-th ket) of the j-th state.

        The coefficients are normalized, i.e. the sum of the absolute values of the coefficients
        in each row is equal to 1.

        """
        return self._cpp.get_coefficients()

    def get_corresponding_ket(self: "Self", state: "StateBase[Any, Any]") -> KetType:
        raise NotImplementedError("Not implemented yet.")

    def get_corresponding_ket_index(self, state: "StateBase[Any, Any]") -> int:
        raise NotImplementedError("Not implemented yet.")

    def get_corresponding_state(self, ket: "KetBase") -> StateType:
        state_cpp = self._cpp.get_corresponding_state(ket._cpp)  # type: ignore [arg-type]
        state_basis = type(self)._from_cpp_object(state_cpp)
        return self._TypeState._from_basis_object(state_basis)

    def get_corresponding_state_index(self, ket: "KetBase") -> int:
        return self._cpp.get_corresponding_state_index(ket._cpp)  # type: ignore [arg-type]

    @abstractmethod
    def get_amplitudes(self, other: Any) -> Any: ...

    @abstractmethod
    def get_overlaps(self, other: Any) -> Any: ...

    @abstractmethod
    def get_matrix_elements(self, other: Any, *args: Any, **kwargs: Any) -> Any: ...
