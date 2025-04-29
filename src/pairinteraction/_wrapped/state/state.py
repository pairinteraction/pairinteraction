# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, Generic, TypeVar

import numpy as np

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.basis import BasisBase
    from pairinteraction._wrapped.ket.ket import KetBase
    from pairinteraction.units import NDArray

KetType = TypeVar("KetType", bound="KetBase", covariant=True)
BasisType = TypeVar("BasisType", bound="BasisBase[Any, Any]", covariant=True)


class StateBase(ABC, Generic[BasisType, KetType]):
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

    _basis: BasisType
    _TypeKet: type[KetType]  # should be ClassVar, but cannot be nested yet

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        raise NotImplementedError(
            "State objects cannot be created directly. "
            "You can use `basis.get_corresponding_state(ket)` or `basis.states[i]` instead."
        )

    @classmethod
    def _from_basis_object(cls: "type[Self]", basis_obj: "BasisBase[KetType, StateBase[Any, Any]]") -> "Self":
        obj = cls.__new__(cls)
        obj._basis = basis_obj  # type: ignore [assignment]
        return obj

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.get_label()})"

    def __str__(self) -> str:
        return self.__repr__().replace("Real", "").replace("Complex", "")

    def get_label(self, max_kets: int = 3) -> str:
        """Label representing the state.

        Returns:
            The label of the ket in the given format.

        """
        coefficients = self.get_coefficients()
        sorted_inds = np.argsort(np.abs(coefficients))[::-1]
        raw = ""
        overlap = 0
        for i, ind in enumerate(sorted_inds, 1):
            raw += f"{coefficients[ind]:.2f} {self.kets[ind].get_label('ket')}"
            overlap += abs(coefficients[ind]) ** 2
            if overlap > 0.9 or i >= max_kets:
                break
            raw += " + "
        if overlap <= 1 - 100 * np.finfo(float).eps:
            raw += " + ... "
        return raw

    @property
    def kets(self) -> list[KetType]:
        """Return a list containing the kets of the basis."""
        return self._basis.kets

    @property
    def number_of_kets(self) -> int:
        """Return the number of kets in the basis."""
        return self._basis.number_of_kets

    def get_coefficients(self) -> "NDArray":
        """Return the coefficients of the state as a 1d-array.

        The coefficients are stored in a numpy.array with shape (number_of_kets,).

        The coefficients are normalized, i.e. the sum of the absolute values of the coefficients is equal to 1.

        """
        return self._basis.get_coefficients().toarray().flatten()

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
