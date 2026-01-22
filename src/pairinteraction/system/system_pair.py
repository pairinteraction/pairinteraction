# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction.basis import BasisPair, BasisPairReal
from pairinteraction.green_tensor import GreenTensorInterpolator
from pairinteraction.system.system_base import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.green_tensor import GreenTensorBase
    from pairinteraction.ket import (
        KetAtom,  # noqa: F401  # needed for sphinx to recognize KetAtomTuple
        KetAtomTuple,
    )
    from pairinteraction.units import (
        ArrayLike,
        NDArray,
        PintArray,  # needed for sphinx to recognize PintArrayLike
        PintArrayLike,
        PintFloat,
    )

logger = logging.getLogger(__name__)


class SystemPair(SystemBase[BasisPair]):
    """System of a pair of atoms.

    Use the given BasisPair object to create the system object.
    You can set the distance (vector) between the atoms afterwards via the corresponding methods.

    Examples:
        >>> import pairinteraction as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> system = pi.SystemAtom(basis).set_magnetic_field([0, 0, 1], unit="G").diagonalize()
        >>> pair_energy = 2 * system.get_corresponding_energy(ket, unit="GHz")
        >>> pair_basis = pi.BasisPair(
        ...     [system, system],
        ...     energy=(pair_energy - 3, pair_energy + 3),
        ...     energy_unit="GHz",
        ... )
        >>> pair_system = pi.SystemPair(pair_basis).set_distance(5, unit="micrometer").set_interaction_order(3)
        >>> print(pair_system)
        SystemPair(BasisPair(|Rb:59,S_1/2,-1/2; Rb:61,S_1/2,-1/2⟩ ... |~Rb:58,F_7/2,7/2; Rb:59,S_1/2,1/2⟩), is_diagonal=False)
        >>> pair_system.diagonalize()
        SystemPair(BasisPair(|Rb:59,S_1/2,-1/2; Rb:61,S_1/2,-1/2⟩ ... |~Rb:58,F_7/2,7/2; Rb:59,S_1/2,1/2⟩), is_diagonal=True)
        >>> eigenenergies = pair_system.get_eigenenergies(unit="GHz")
        >>> print(f"{eigenenergies[0] - pair_energy:.5f}")
        -2.18394

    """  # noqa: E501

    _cpp: _backend.SystemPairComplex
    _cpp_type = _backend.SystemPairComplex
    _basis_class = BasisPair

    def __init__(self, basis: BasisPair) -> None:
        """Create a system object for a pair of atoms.

        Args:
            basis: The :class:`pairinteraction.BasisPair` object that describes the basis of the system.

        """
        self._cpp = self._cpp_type(basis._cpp)
        self._basis = basis
        self._distance_vector_au = np.array([0, 0, np.inf])
        self._interaction_order = 3

    def get_eigenbasis(self) -> BasisPair:
        """Get the eigenbasis of the system.

        This method retrieves the eigenbasis of the system, which is the basis in which the Hamiltonian is diagonal.

        Returns:
            A BasisPair object representing the eigenbasis of the system.

        """
        eigenbasis = super().get_eigenbasis()
        eigenbasis.system_atoms = self.basis.system_atoms
        return eigenbasis

    def set_interaction_order(self: Self, order: int) -> Self:
        """Set the interaction order of the pair system.

        Args:
            order: The interaction order to set for the pair system.
                The order must be 3, 4, or 5.

        """
        self._interaction_order = order
        self._cpp.set_interaction_order(order)
        return self

    def set_distance(
        self: Self,
        distance: float | PintFloat,
        angle_degree: float = 0,
        unit: str | None = None,
    ) -> Self:
        """Set the distance between the atoms using the specified distance and angle.

        Args:
            distance: The distance to set between the atoms in the given unit.
            angle_degree: The angle between the distance vector and the z-axis in degrees.
                90 degrees corresponds to the x-axis.
                Defaults to 0, which corresponds to the z-axis.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        distance_vector = [np.sin(np.deg2rad(angle_degree)) * distance, 0, np.cos(np.deg2rad(angle_degree)) * distance]
        return self.set_distance_vector(distance_vector, unit)

    def set_distance_vector(
        self: Self,
        distance: ArrayLike | PintArrayLike,
        unit: str | None = None,
    ) -> Self:
        """Set the distance vector between the atoms.

        Args:
            distance: The distance vector to set between the atoms in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        distance_au = [QuantityScalar.convert_user_to_au(v, unit, "distance") for v in distance]
        self._cpp.set_distance_vector(distance_au)
        self._distance_vector_au = np.array(distance_au)
        return self

    @overload
    def get_distance_vector(self, unit: None = None) -> PintArray: ...

    @overload
    def get_distance_vector(self, unit: str) -> NDArray: ...

    def get_distance_vector(self, unit: str | None = None) -> NDArray | PintArray:
        return np.array([QuantityScalar.convert_au_to_user(d, "distance", unit) for d in self._distance_vector_au])  # type: ignore [return-value]

    @overload
    def get_distance(self, unit: None = None) -> PintFloat: ...

    @overload
    def get_distance(self, unit: str) -> float: ...

    def get_distance(self, unit: str | None = None) -> float | PintFloat:
        distance = np.linalg.norm(self._distance_vector_au)
        return QuantityScalar.convert_au_to_user(float(distance), "distance", unit)

    def set_green_tensor_interpolator(self, green_tensor_interpolator: GreenTensorInterpolator) -> Self:
        """Set the Green tensor interpolator for the pair system.

        Args:
            green_tensor_interpolator: The Green tensor interpolator to set for the system.

        """
        self._cpp.set_green_tensor_interpolator(green_tensor_interpolator._cpp)
        return self

    def set_green_tensor(self, green_tensor: GreenTensorBase, omega_steps: int) -> Self:
        """Set the Green tensor for the pair system.

        Args:
            green_tensor: The Green tensor to set for the system.
            omega_steps: The number of omega steps to use for non-constant Green tensors.
                If 1, a constant Green tensor (with omega=0) is assumed.

        """
        if green_tensor.pos1_au is None or green_tensor.pos2_au is None:
            raise ValueError("The positions of the atoms in the Green tensor must be set before using it.")
        self._distance_vector_au = green_tensor.pos1_au - green_tensor.pos2_au

        if omega_steps <= 0:
            raise ValueError("omega_steps must be a positive integer.")

        if omega_steps == 1:
            gt = green_tensor.get_dipole_dipole(omega=1e-10, omega_unit="Hz")
            gti = GreenTensorInterpolator()
            gti.set_from_cartesian(1, 1, gt.magnitude, tensor_unit=str(gt.units))
            self.set_green_tensor_interpolator(gti)
            return self

        # TODO optimize how to choose omegas, for now we just use linear spacing
        energies_au = [ket.get_energy("hartree") for ket in self.basis.kets]
        omega_max = max(energies_au) - min(energies_au)

        gti = green_tensor.get_green_tensor_interpolator(1e-10, omega_max, omega_steps, omega_unit="hartree")
        self.set_green_tensor_interpolator(gti)
        return self

    @overload
    def get_corresponding_energy(self: Self, ket_tuple: KetAtomTuple, unit: None = None) -> PintFloat: ...

    @overload
    def get_corresponding_energy(self: Self, ket_tuple: KetAtomTuple, unit: str) -> float: ...

    def get_corresponding_energy(self: Self, ket_tuple: KetAtomTuple, unit: str | None = None) -> float | PintFloat:
        overlaps = self.get_eigenbasis().get_overlaps(ket_tuple)
        idx = np.argmax(overlaps)
        if overlaps[idx] <= 0.5:
            logger.warning(
                "The provided ket states does not correspond to an eigenstate of the system in a unique way."
            )
        return self.get_eigenenergies(unit=unit)[idx]  # type: ignore [index,no-any-return] # PintArray does not know it can be indexed


class SystemPairReal(SystemPair):
    _cpp: _backend.SystemPairReal  # type: ignore [assignment]
    _cpp_type = _backend.SystemPairReal  # type: ignore [assignment]
    _basis_class = BasisPairReal
