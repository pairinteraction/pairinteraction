# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const
from typing_extensions import override

from pairinteraction.green_tensor.dynamic_green_tensor import (
    dynamic_green_tensor_homogeneous,
    dynamic_green_tensor_scattered,
)
from pairinteraction.green_tensor.green_tensor_base import (
    GreenTensorBase,
    _get_lab_to_local_rotation_matrix,
    _rotate_tensor_to_lab,
    _rotate_vector_to_local,
    evaluate_relative_permittivity,
)
from pairinteraction.units import QuantityScalar, ureg

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.green_tensor.green_tensor_base import PermittivityLike
    from pairinteraction.units import ArrayLike, NDArray, PintArrayLike


class GreenTensorSurface(GreenTensorBase):
    """Green tensor for two atoms near a single infinite surface.

    Examples:
        >>> from pairinteraction.green_tensor import GreenTensorSurface
        >>> gt = GreenTensorSurface([0, 0, 0], [10, 0, 0], point_on_plane=[0, 0, -5], unit="micrometer")
        >>> transition_energy = 2  # h * GHz
        >>> gt_dipole_dipole = gt.get(1, 1, transition_energy, "planck_constant * GHz")
        >>> print(f"{gt_dipole_dipole[0, 0]:.2f}")
        -4.37 / bohr

    """

    def __init__(
        self,
        pos1: ArrayLike | PintArrayLike,
        pos2: ArrayLike | PintArrayLike,
        point_on_plane: ArrayLike | PintArrayLike,
        normal: ArrayLike | PintArrayLike = (0, 0, 1),
        unit: str | None = None,
        static_limit: bool = True,
        interaction_order: int = 3,
        *,
        without_vacuum_contribution: bool = False,
    ) -> None:
        """Create a Green tensor for two atoms near a single infinite surface.

        The surface is an infinite plane specified by a point on the plane and its normal vector.
        If not specified otherwise (see `set_relative_permittivities`), the surface is treated as a perfect mirror.


        Args:
            pos1: Position of the first atom in the given unit.
            pos2: Position of the second atom in the given unit.
            point_on_plane: A point on the surface plane in the given unit.
            normal: The surface normal vector. Defaults to `[0, 0, 1]`.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.
            static_limit: If True, the static limit is used.
                Default True.
            interaction_order: The order of interaction, e.g., 3 for dipole-dipole.
                Defaults to 3.
            without_vacuum_contribution: If True, return only the scattered contribution and
                omit the homogeneous vacuum term. Defaults to False.

        """
        super().__init__(
            pos1, pos2, unit, static_limit, interaction_order, without_vacuum_contribution=without_vacuum_contribution
        )
        self.point_on_plane_au = np.array(
            [QuantityScalar.convert_user_to_au(v, unit, "distance") for v in point_on_plane]
        )
        self.normal_au = np.array(normal, dtype=float)
        self._lab_to_local_rotation = _get_lab_to_local_rotation_matrix(self.normal_au)
        self.surface_epsilon: PermittivityLike = 1e9  # Almost perfect mirror # TODO make utils be able to handle inf

    def set_relative_permittivities(self, epsilon: PermittivityLike, surface_epsilon: PermittivityLike) -> Self:
        """Set the relative permittivities of the system.

        Args:
            epsilon: The relative permittivity (dimensionless) of the medium inside the cavity.
            surface_epsilon: The relative permittivity (dimensionless) of the surface.


        """
        if self.without_vacuum_contribution and epsilon != 1.0:
            raise ValueError(
                "If the Green tensor is provided without the vacuum contribution, "
                "the relative permittivity of the medium needs to be set to 1."
            )

        self.epsilon = epsilon
        self.surface_epsilon = surface_epsilon
        return self

    @override
    def _get_scaled_au(self, kappa1: int, kappa2: int, transition_energy_au: float) -> NDArray:
        if kappa1 == 1 and kappa2 == 1:
            return self._get_scaled_dipole_dipole_au(transition_energy_au)
        raise NotImplementedError("Only dipole-dipole Green tensors are currently implemented.")

    def _get_scaled_dipole_dipole_au(self, transition_energy_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for a single surface in atomic units.

        Args:
            transition_energy_au: The transition energy in atomic units at which to evaluate the Green tensor.

        Returns:
            The dipole dipole Green tensor in cartesian coordinates as a 3x3 array in atomic units (i.e. 1/bohr).

        """
        au_to_meter: float = ureg.Quantity(1, "atomic_unit_of_length").to("meter").magnitude
        pos1_local_m = _rotate_vector_to_local(np.array(self.pos1_au) * au_to_meter, self._lab_to_local_rotation)
        pos2_local_m = _rotate_vector_to_local(np.array(self.pos2_au) * au_to_meter, self._lab_to_local_rotation)
        point_on_plane_local_m = _rotate_vector_to_local(
            self.point_on_plane_au * au_to_meter, self._lab_to_local_rotation
        )
        epsilon = evaluate_relative_permittivity(self.epsilon, transition_energy_au, "hartree")

        omega_hz = ureg.Quantity(transition_energy_au, "hartree").to("hbar Hz", "spectroscopy").magnitude

        # Assume two surfaces, where the further apart atom is located in the center
        # but the second surface has the same permittivity as the inbetween medium
        z1_m = point_on_plane_local_m[2]
        height = 2 * max(abs(pos1_local_m[2] - z1_m), abs(pos2_local_m[2] - z1_m))
        if pos1_local_m[2] < z1_m and pos2_local_m[2] < z1_m:
            z2_m = z1_m - height
        elif pos1_local_m[2] > z1_m and pos2_local_m[2] > z1_m:
            z2_m = z1_m + height
        else:
            raise ValueError("Both atoms must be located either above or below the surface.")

        epsilon1 = evaluate_relative_permittivity(self.surface_epsilon, transition_energy_au, "hartree")
        epsilon2 = epsilon

        # unit: # m^(-3) [hbar]^(-1) [epsilon_0]^(-1)
        gt = dynamic_green_tensor_scattered(
            pos1_local_m, pos2_local_m, z1_m, z2_m, omega_hz, epsilon, epsilon1, epsilon2, only_real_part=True
        )
        if not self.without_vacuum_contribution:
            gt += dynamic_green_tensor_homogeneous(pos1_local_m, pos2_local_m, omega_hz, epsilon, only_real_part=True)
        gt = _rotate_tensor_to_lab(gt, self._lab_to_local_rotation)
        to_au = au_to_meter ** (-3) * ((4 * np.pi) ** (-1)) / (const.epsilon_0 * const.hbar)
        # hbar * epsilon_0 = (4*np.pi)**(-1) in atomic units
        return np.real(gt) / to_au
