# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const

from pairinteraction.green_tensor import utils
from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase, get_electric_permitivity
from pairinteraction.units import QuantityScalar, ureg

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.green_tensor.green_tensor_base import Permitivity
    from pairinteraction.units import ArrayLike, NDArray, PintArrayLike, PintFloat


class GreenTensorSurface(GreenTensorBase):
    def __init__(
        self,
        pos1: ArrayLike | PintArrayLike,
        pos2: ArrayLike | PintArrayLike,
        z: float | PintFloat,
        unit: str | None = None,
        static_limit: bool = False,
        interaction_order: int = 3,
    ) -> None:
        """Create a Green tensor for two atoms near a single infinite surface.

        The surface is assumed to be infinite in the x-y plane.

        Args:
            pos1: Position of the first atom in the given unit.
            pos2: Position of the second atom in the given unit.
            z: The z-position of the surface in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.
            static_limit: If True, the static limit is used.
                Default False.
            interaction_order: The order of interaction, e.g., 3 for dipole-dipole.
                Defaults to 3.

        """
        super().__init__(pos1, pos2, unit, static_limit, interaction_order)
        self.surface_z_au = QuantityScalar.convert_user_to_au(z, unit, "distance")

    def set_relative_permittivities(self, epsilon: Permitivity, surface_epsilon: Permitivity) -> Self:
        """Set the relative permittivities of the system.

        Args:
            epsilon: The relative permittivity (dimensionless) of the medium inside the cavity.
            surface_epsilon: The relative permittivity (dimensionless) of the surface.


        """
        self.epsilon = epsilon
        self.surface_epsilon = surface_epsilon
        return self

    def _get_scaled_dipole_dipole_au(self, omega_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for a single surface in atomic units.

        Args:
            omega_au: The angular frequency in atomic units at which to evaluate the Green tensor.

        Returns:
            The dipole dipole Green tensor in cartesian coordinates as a 3x3 array in atomic units (i.e. 1/bohr).

        """
        if self.pos1_au is None or self.pos2_au is None:
            raise RuntimeError("Atom positions have to be set before calculating the Green tensor.")

        au_to_meter: float = ureg.Quantity(1, "atomic_unit_of_length").to("meter").magnitude
        pos1_m = np.array(self.pos1_au) * au_to_meter
        pos2_m = np.array(self.pos2_au) * au_to_meter
        epsilon = get_electric_permitivity(self.epsilon, omega_au, "hartree")

        omega = ureg.Quantity(omega_au, "hartree").to("Hz", "spectroscopy")
        omega_hz = omega.magnitude

        surface_z = self.surface_z_au * au_to_meter

        # Assume two surfaces, but the second surface has the same permittivity as the inbetween medium
        height = 2 * max(abs(surface_z - pos1_m[2]), abs(surface_z - pos2_m[2]))
        if pos1_m[2] < surface_z and pos2_m[2] < surface_z:
            surface2_z = surface_z - height
        elif pos1_m[2] > surface_z and pos2_m[2] > surface_z:
            surface2_z = surface_z + height
        else:
            raise ValueError("Both atoms must be located either above or below the surface.")

        pos1_shifted_m = pos1_m - np.array([0, 0, min(surface_z, surface2_z)])
        pos2_shifted_m = pos2_m - np.array([0, 0, min(surface_z, surface2_z)])

        if surface2_z < surface_z:
            epsilon_top = get_electric_permitivity(self.surface_epsilon, omega_au, "hartree")
            epsilon_bottom = epsilon
        else:
            epsilon_top = epsilon
            epsilon_bottom = get_electric_permitivity(self.surface_epsilon, omega_au, "hartree")

        # unit: # m^(-3) [hbar]^(-1) [epsilon_0]^(-1)
        gt = utils.green_tensor_total(
            pos1_shifted_m, pos2_shifted_m, omega_hz, epsilon, epsilon_top, epsilon_bottom, height, only_real_part=True
        )
        to_au = au_to_meter ** (-3) * ((4 * np.pi) ** (-1)) / (const.epsilon_0 * const.hbar)
        # hbar * epsilon_0 = (4*np.pi)**(-1) in atomc units
        return np.real(gt) / to_au
