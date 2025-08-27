# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from abc import ABC
from typing import TYPE_CHECKING, Any, ClassVar, Generic, Optional, TypeVar, Union, overload

import numpy as np
from typing_extensions import deprecated

from pairinteraction import _backend
from pairinteraction._wrapped.diagonalize.diagonalize import diagonalize
from pairinteraction._wrapped.diagonalize.diagonalizer import Diagonalizer
from pairinteraction._wrapped.enums import FloatType
from pairinteraction.units import UnitConverterArray, UnitConverterSparse

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.basis import BasisBase
    from pairinteraction.units import NDArray, PintArray, PintFloat, PintSparse

    Quantity = TypeVar("Quantity", float, "PintFloat")


BasisType = TypeVar("BasisType", bound="BasisBase[Any, Any]", covariant=True)
UnionCPPSystem = Union[
    _backend.SystemAtomReal, _backend.SystemAtomComplex, _backend.SystemPairReal, _backend.SystemPairComplex
]
UnionTypeCPPSystem = Union[
    type[_backend.SystemAtomReal],
    type[_backend.SystemAtomComplex],
    type[_backend.SystemPairReal],
    type[_backend.SystemPairComplex],
]


class SystemBase(ABC, Generic[BasisType]):
    """Base class for all System objects.

    The system objects are meant to represent the full physical system, including the basis
    as well as possible fields and interactions which form the Hamiltonian.

    Typically, the system objects are simply constructed from a basis object.
    The fields and interactions can then be set afterwards.

    All system objects share a few common attributes and methods, that are defined in this base class, e.g.:
        - the basis of the system,
        - the matrix stored as scipy sparse matrix, which represents the Hamiltonian,
        - ...

    """

    _cpp: UnionCPPSystem
    _cpp_type: ClassVar[UnionTypeCPPSystem]
    _basis: BasisType
    _TypeBasis: type[BasisType]  # should be ClassVar, but cannot be nested yet

    def __init__(self, basis: BasisType) -> None:
        self._cpp = self._cpp_type(basis._cpp)  # type: ignore [arg-type]
        self._basis = basis

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.basis!r}, is_diagonal={self.is_diagonal})"

    def __str__(self) -> str:
        return self.__repr__().replace("Real", "").replace("Complex", "")

    def _update_basis(self) -> None:
        self._basis = self._TypeBasis._from_cpp_object(self._cpp.get_basis())

    @overload
    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "eigen",
        float_type: FloatType = "float64",
        rtol: float = 1e-6,
        sort_by_energy: bool = True,
        energy_range: tuple[Union["Quantity", None], Union["Quantity", None]] = (None, None),
        energy_range_unit: Optional[str] = None,
        m0: Optional[int] = None,
    ) -> "Self": ...

    @overload
    @deprecated("Use energy_range_unit=... instead of energy_unit=...")
    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "eigen",
        float_type: FloatType = "float64",
        rtol: float = 1e-6,
        sort_by_energy: bool = True,
        energy_range: tuple[Union["Quantity", None], Union["Quantity", None]] = (None, None),
        *,
        energy_unit: Optional[str],
        m0: Optional[int] = None,
    ) -> "Self": ...

    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "eigen",
        float_type: FloatType = "float64",
        rtol: float = 1e-6,
        sort_by_energy: bool = True,
        energy_range: tuple[Union["Quantity", None], Union["Quantity", None]] = (None, None),
        energy_range_unit: Optional[str] = None,
        m0: Optional[int] = None,
        *,
        energy_unit: Optional[str] = None,
    ) -> "Self":
        """Diagonalize the Hamiltonian and update the basis to the eigenbasis.

        This method diagonalizes the Hamiltonian of the system and updates the basis of the system accordingly.

        Args:
            systems: A list of `SystemAtom` or `SystemPair` objects.
            diagonalizer: The diagonalizer method to use. Defaults to "eigen".
            float_type: The floating point precision to use for the diagonalization. Defaults to "float64".
            rtol: The relative tolerance allowed for eigenenergies. The error in eigenenergies is bounded
                by rtol * ||H||, where ||H|| is the norm of the Hamiltonian matrix. Defaults to 1e-6.
            sort_by_energy: Whether to sort the resulting basis by energy. Defaults to True.
            energy_range: A tuple specifying an energy range, in which the eigenenergies should be calculated.
                Specifying a range can speed up the diagonalization process (depending on the diagonalizer method).
                The accuracy of the eigenenergies is not affected by this, but not all eigenenergies will be calculated.
                Defaults to (None, None), i.e. calculate all eigenenergies.
            energy_range_unit: The unit in which the energy_range is given. Defaults to None assumes pint objects.
            m0: The search subspace size for the FEAST diagonalizer. Defaults to None.
            energy_unit: Deprecated, use energy_range_unit instead.

        Returns:
            Self: The updated instance of the system.

        """
        diagonalize(
            [self],
            diagonalizer,
            float_type,
            rtol,
            sort_by_energy,
            energy_range,
            energy_range_unit,
            m0,
            energy_unit=energy_unit,
        )  # type: ignore [misc,call-overload]
        return self

    @property
    def is_diagonal(self) -> bool:
        return self._cpp.is_diagonal()

    @property
    def basis(self) -> BasisType:
        return self._basis

    @property
    def matrix(self) -> "csr_matrix":
        return self._cpp.get_matrix()

    def get_eigenbasis(self) -> BasisType:
        cpp_eigenbasis = self._cpp.get_eigenbasis()
        return self._TypeBasis._from_cpp_object(cpp_eigenbasis)

    @overload
    def get_eigenenergies(self, unit: None = None) -> "PintArray": ...

    @overload
    def get_eigenenergies(self, unit: str) -> "NDArray": ...

    def get_eigenenergies(self, unit: Optional[str] = None) -> Union["NDArray", "PintArray"]:
        eigenenergies_au: NDArray = np.array(self._cpp.get_eigenenergies())
        return UnitConverterArray.au_to_user(eigenenergies_au, "energy", unit)

    @overload
    def get_hamiltonian(self, unit: None = None) -> "PintSparse": ...  # type: ignore [type-var] # see PintSparse

    @overload
    def get_hamiltonian(self, unit: str) -> "csr_matrix": ...

    def get_hamiltonian(self, unit: Optional[str] = None) -> Union["csr_matrix", "PintSparse"]:
        hamiltonian_au = self._cpp.get_matrix()
        return UnitConverterSparse.au_to_user(hamiltonian_au, "energy", unit)
