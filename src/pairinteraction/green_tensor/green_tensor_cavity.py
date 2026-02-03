# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np
import scipy.constants as const

from pairinteraction.green_tensor import utils
from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase, get_electric_permitivity
from pairinteraction.units import QuantityScalar, ureg

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.units import NDArray, PintFloat


class GreenTensorCavity(GreenTensorBase):
    def set_surface_positions(self, z1: float, z2: float, unit: str | None = None) -> Self:
        """Set the positions of the two surfaces of the cavity along the z-axis.

        The two surfaces of the cavity are assumed to be infinite in the x-y plane.

        Args:
            z1: The z-position of the first surface in the given unit.
            z2: The z-position of the second surface in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        self.surface1_z_au = QuantityScalar.convert_user_to_au(z1, unit, "distance")
        self.surface2_z_au = QuantityScalar.convert_user_to_au(z2, unit, "distance")
        return self

    def set_electric_permitivity_surfaces(
        self, epsilon1: complex | Callable[[PintFloat], complex], epsilon2: complex | Callable[[PintFloat], complex]
    ) -> Self:
        """Set the electric permittivity for the surfaces.

        Args:
            epsilon1: The electric permittivity (dimensionless) of the first surface.
            epsilon2: The electric permittivity (dimensionless) of the second surface.

        """
        self.surface1_epsilon = epsilon1
        self.surface2_epsilon = epsilon2
        return self

    def _get_scaled_dipole_dipole_au(self, omega_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for a cavity in atomic units.

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

        surface1_z_m = self.surface1_z_au * au_to_meter
        surface2_z_m = self.surface2_z_au * au_to_meter

        height = abs(surface1_z_m - surface2_z_m)
        pos1_shifted_m = pos1_m - np.array([0, 0, min(surface1_z_m, surface2_z_m)])
        pos2_shifted_m = pos2_m - np.array([0, 0, min(surface1_z_m, surface2_z_m)])

        if not 0 < pos1_shifted_m[2] < height or not 0 < pos2_shifted_m[2] < height:
            raise ValueError("Both atoms must be located inside the cavity between the two surfaces.")

        if surface2_z_m < surface1_z_m:
            epsilon_top = get_electric_permitivity(self.surface1_epsilon, omega_au, "hartree")
            epsilon_bottom = get_electric_permitivity(self.surface2_epsilon, omega_au, "hartree")
        else:
            epsilon_top = get_electric_permitivity(self.surface2_epsilon, omega_au, "hartree")
            epsilon_bottom = get_electric_permitivity(self.surface1_epsilon, omega_au, "hartree")

        # unit: # m^(-3) [hbar]^(-1) [epsilon_0]^(-1)
        gt = utils.green_tensor_total(
            pos1_shifted_m, pos2_shifted_m, omega_hz, epsilon, epsilon_top, epsilon_bottom, height, only_real_part=True
        )
        to_au = au_to_meter**(-3) * ((4 * np.pi) ** (-1)) / (const.epsilon_0 * const.hbar)
        # hbar * epsilon_0 = (4*np.pi)**(-1) in atomc units
        return np.real(gt) / to_au
