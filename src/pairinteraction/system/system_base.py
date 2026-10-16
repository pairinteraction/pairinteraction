# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING, Any, Generic, TypeAlias, TypeVar, overload

import numpy as np

from pairinteraction.basis.basis_pair import BasisPair
from pairinteraction.diagonalization import diagonalize
from pairinteraction.state.state_base import get_index_with_largest_overlap
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction import _backend
    from pairinteraction.basis import BasisBase
    from pairinteraction.diagonalization import Diagonalizer
    from pairinteraction.enums import FloatType
    from pairinteraction.ket import KetBase
    from pairinteraction.ket.ket_pair import KetPairLike
    from pairinteraction.units import NDArray, PintArray, PintFloat, PintSparse

    Quantity = TypeVar("Quantity", float, "PintFloat")


BasisType = TypeVar("BasisType", bound="BasisBase[Any, Any]")
UnionCPPSystem: TypeAlias = "_backend.SystemAtomComplex | _backend.SystemPairComplex"


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
    _basis: BasisType
    _basis_class: type[BasisType]  # should be ClassVar, but cannot be nested yet

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.basis!r}, is_diagonal={self.is_diagonal})"

    def __str__(self) -> str:
        return self.__repr__()

    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "eigen",
        float_type: FloatType = "float64",
        rtol: float = 1e-6,
        sort_by_energy: bool = True,
        energy_range: tuple[Quantity | None, Quantity | None] = (None, None),
        energy_range_unit: str | None = None,
        m0: int | None = None,
    ) -> Self:
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

        Returns:
            Self: The updated instance of the system.

        """
        diagonalize([self], diagonalizer, float_type, rtol, sort_by_energy, energy_range, energy_range_unit, m0)
        return self

    @property
    def is_diagonal(self) -> bool:
        return self._cpp.is_diagonal()

    @property
    def basis(self) -> BasisType:
        """The basis object of the system."""
        cpp_basis = self._cpp.get_basis()
        if self._basis._cpp != cpp_basis:
            kwargs = {"system_atoms": self._basis.system_atoms} if isinstance(self._basis, BasisPair) else {}
            self._basis = self._basis_class._from_cpp_object(cpp_basis, **kwargs)
        return self._basis

    @property
    def matrix(self) -> csr_matrix:
        return self._cpp.get_matrix()

    def get_eigenbasis(self) -> BasisType:
        _eigenbasis: BasisType | None = getattr(self, "_eigenbasis", None)
        cpp_eigenbasis = self._cpp.get_eigenbasis()
        if _eigenbasis is None or _eigenbasis._cpp != cpp_eigenbasis:
            kwargs = {"system_atoms": self._basis.system_atoms} if isinstance(self._basis, BasisPair) else {}
            self._eigenbasis = self._basis_class._from_cpp_object(cpp_eigenbasis, **kwargs)
        return self._eigenbasis

    @overload
    def get_eigenenergies(self, unit: None = None) -> PintArray: ...

    @overload
    def get_eigenenergies(self, unit: str) -> NDArray: ...

    def get_eigenenergies(self, unit: str | None = None) -> NDArray | PintArray:
        eigenenergies_au: NDArray = np.array(self._cpp.get_eigenenergies())
        return QuantityArray.convert_au_to_user(eigenenergies_au, "energy", unit)

    @overload
    def get_corresponding_energy(self, ket: KetBase | KetPairLike, unit: None = None) -> PintFloat: ...

    @overload
    def get_corresponding_energy(self, ket: KetBase | KetPairLike, unit: str) -> float: ...

    def get_corresponding_energy(self, ket: KetBase | KetPairLike, unit: str | None = None) -> float | PintFloat:
        """Return the energy of the eigenstate with the largest overlap with the given ket.

        If the ket is distributed over several degenerate eigenstates (e.g. symmetrized pair states),
        the corresponding energy is still unique, so no warning is logged in this case.
        """
        overlaps = self.get_eigenbasis().get_overlaps(ket)
        eigenenergies_au: NDArray = np.array(self._cpp.get_eigenenergies())

        state_idx = int(np.argmax(overlaps))
        tolerance = 100 * np.finfo(float).eps
        energy_au = eigenenergies_au[state_idx]
        is_degenerate = np.abs(eigenenergies_au - energy_au) <= tolerance * abs(energy_au)
        if np.sum(overlaps[is_degenerate]) <= 0.5 + tolerance:
            # The ket is not mainly in one (degenerate) eigenspace, call only to log the warning
            get_index_with_largest_overlap(overlaps, err_msg="energy of an eigenstate for the given ket")

        return QuantityScalar.convert_au_to_user(energy_au, "energy", unit)

    @overload
    def get_hamiltonian(self, unit: None = None) -> PintSparse: ...

    @overload
    def get_hamiltonian(self, unit: str) -> csr_matrix: ...

    def get_hamiltonian(self, unit: str | None = None) -> csr_matrix | PintSparse:
        hamiltonian_au = self._cpp.get_matrix()
        hamiltonian_au.data = np.real_if_close(hamiltonian_au.data)
        return QuantitySparse.convert_au_to_user(hamiltonian_au, "energy", unit)
