# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Generic, TypeAlias, TypeVar

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction import _backend
    from pairinteraction.ket import KetBase
    from pairinteraction.state import StateBase

KetType = TypeVar("KetType", bound="KetBase")
StateType = TypeVar("StateType", bound="StateBase[Any]")
UnionCPPBasis: TypeAlias = "_backend.BasisAtomComplex | _backend.BasisPairComplex"


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
    _ket_class: type[KetType]  # should be ClassVar, but cannot be nested yet
    _state_class: type[StateType]  # should be ClassVar, but cannot be nested yet

    def _post_init(self) -> None:
        if self.number_of_kets == 0:
            raise ValueError("Cannot create a basis with zero kets.")
        self._kets_cache: dict[int, KetType] = {}

    @classmethod
    def _from_cpp_object(cls: type[Self], cpp_obj: UnionCPPBasis, *args: Any, **kwargs: Any) -> Self:
        assert len(kwargs) == len(args) == 0, "No additional arguments expected."
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        obj._post_init()
        return obj

    def __repr__(self) -> str:
        args = f"{self.get_ket(0)} ... {self.get_ket(self.number_of_kets - 1)}"
        return f"{type(self).__name__}({args})"

    def __str__(self) -> str:
        return self.__repr__()

    @property
    def kets(self) -> list[KetType]:
        """Return a list containing the kets of the basis."""
        return [self.get_ket(i) for i in range(self.number_of_kets)]

    def get_ket(self, index: int) -> KetType:
        """Return the ket at the given index."""
        if index not in self._kets_cache:
            if index < 0 or index >= self.number_of_kets:
                raise IndexError(f"Ket index {index} out of range (number of kets: {self.number_of_kets}).")
            ket_cpp = self._cpp.get_ket(index)
            self._kets_cache[index] = self._ket_class._from_cpp_object(ket_cpp)
        return self._kets_cache[index]

    @property  # not cached_property since State objects are mutable
    def states(self) -> list[StateType]:
        """Return a list containing the states of the basis."""
        return [self.get_state(i) for i in range(self.number_of_states)]

    def get_state(self, index: int) -> StateType:
        """Return the state at the given index."""
        if index < 0 or index >= self.number_of_states:
            raise IndexError(f"State index {index} out of range (number of states: {self.number_of_states}).")
        state_cpp = self._cpp.get_state(index)
        return self._state_class._from_cpp_object(state_cpp)

    @property
    def number_of_kets(self) -> int:
        """Return the number of kets in the basis."""
        return self._cpp.get_number_of_kets()

    @property
    def number_of_states(self) -> int:
        """Return the number of states in the basis."""
        return self._cpp.get_number_of_states()

    def get_coefficients(self) -> csr_matrix:
        """Return the coefficients of the basis as a sparse matrix.

        The coefficients are stored in a sparse matrix with shape (number_of_kets, number_of_states),
        where the first index correspond to the kets and the second index correspond to the states.
        For example `basis.get_coefficients()[i, j]` is the i-th coefficient
        (i.e. the coefficient corresponding to the i-th ket) of the j-th state.

        The coefficients are normalized, i.e. the sum of the absolute values of the coefficients
        in each row is equal to 1.

        """
        return self._cpp.get_coefficients()

    def get_corresponding_ket(self: Self, state: StateType) -> KetType:
        ket_index = self.get_corresponding_ket_index(state)
        return self.get_ket(ket_index)

    def get_corresponding_ket_index(self, state: StateType) -> int:
        raise NotImplementedError("Not implemented yet.")

    def get_corresponding_state(self, ket: KetBase) -> StateType:
        state_index = self.get_corresponding_state_index(ket)
        return self.get_state(state_index)

    def get_corresponding_state_index(self, ket: KetBase) -> int:
        return self._cpp.get_corresponding_state_index(ket._cpp)  # type: ignore [arg-type]

    def canonicalized(self: Self) -> Self:
        """Return the canonical basis with identity coefficients."""
        return type(self)._from_cpp_object(self._cpp.canonicalized())

    @abstractmethod
    def get_amplitudes(self, other: Any) -> Any: ...

    @abstractmethod
    def get_overlaps(self, other: Any) -> Any: ...

    @abstractmethod
    def get_matrix_elements(self, other: Any, *args: Any, **kwargs: Any) -> Any: ...
