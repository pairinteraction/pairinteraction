# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Generic, TypeVar, Union

import numpy as np

from pairinteraction import _backend

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.ket import KetBase
    from pairinteraction.units import NDArray

KetType = TypeVar("KetType", bound="KetBase")
UnionCPPBasis = Union[_backend.BasisAtomComplex, _backend.BasisPairComplex]


class StateBase(ABC, Generic[KetType]):
    """Base class for all State objects.

    The state objects are meant to represent a set of kets, that span a Hilbert space
    and store a coefficient vector, that describes the state in terms of the kets.

    Basically just a wrapper around the python Basis classes, but with more convenience functions specific to a state.

    Note, in the cpp code we dont have an explicit state class, but use a Basis object, which has only one state

    All state objects share a few common attributes and methods, that are defined in this base class, e.g.:
        - the number of kets,
        - the kets of the state,
        - the coefficient vector as 1d-array,
        - ...
    """

    _cpp: UnionCPPBasis
    _ket_class: type[KetType]  # should be ClassVar, but cannot be nested yet

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        raise NotImplementedError(
            "State objects cannot be created directly. "
            "You can use `basis.get_corresponding_state(ket)` or `basis.states[i]` instead."
        )

    @classmethod
    def _from_cpp_object(cls: type[Self], cpp_obj: UnionCPPBasis) -> Self:
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.get_label()})"

    def __str__(self) -> str:
        return self.__repr__()

    def get_label(self, max_kets: int = 3) -> str:
        """Label representing the state.

        Args:
            max_kets: Maximum number of kets to include in the label.

        Returns:
            The label of the ket in the given format.

        """
        coefficients = self.get_coefficients()
        sorted_inds = np.argsort(np.abs(coefficients))[::-1]
        norm_squared = self.norm**2
        label = ""
        overlap = 0
        for i, ind in enumerate(sorted_inds, 1):
            label += f"{np.real_if_close(coefficients[ind]):.2f} {self.kets[ind].get_label('ket')}"
            overlap += abs(coefficients[ind]) ** 2
            if overlap > (0.95 * norm_squared) or i >= max_kets:
                break
            label += " + "
        if overlap <= norm_squared - 100 * np.finfo(float).eps:
            label += " + ... "
        return label

    @property
    def kets(self) -> list[KetType]:
        """Return a list containing the kets of the basis."""
        return [self._ket_class._from_cpp_object(ket_cpp) for ket_cpp in self._cpp.get_kets()]

    @property
    def number_of_kets(self) -> int:
        """Return the number of kets in the basis."""
        return self._cpp.get_number_of_kets()

    @property
    def norm(self) -> np.floating:
        """Return the norm of the state."""
        return np.linalg.norm(self.get_coefficients())

    def get_coefficients(self) -> NDArray:
        """Return the coefficients of the state as a 1d-array.

        The coefficients are stored in a numpy.array with shape (number_of_kets,).

        The coefficients are normalized, i.e. the sum of the absolute values of the coefficients is equal to 1.

        """
        return self._cpp.get_coefficients().toarray().flatten()

    def get_corresponding_ket(self) -> KetType:
        """Return the ket corresponding to the state (i.e. the ket with the maximal overlap)."""
        return self.kets[self.get_corresponding_ket_index()]

    def get_corresponding_ket_index(self) -> int:
        """Return the index of the ket corresponding to the state (i.e. the ket with the maximal overlap)."""
        coeffs = np.abs(self.get_coefficients())
        ket_idx = np.argmax(coeffs)
        if coeffs[ket_idx] ** 2 < 0.5 + 100 * np.finfo(float).eps:
            raise ValueError("The maximal overlap is <= 0.5, i.e. the state has no uniquely corresponding ket.")
        return int(ket_idx)

    @abstractmethod
    def get_amplitude(self, other: Any) -> Any: ...

    @abstractmethod
    def get_overlap(self, other: Any) -> Any: ...

    @abstractmethod
    def get_matrix_element(self, other: Any, *args: Any, **kwargs: Any) -> Any: ...
