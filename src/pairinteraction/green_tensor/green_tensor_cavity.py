# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from green_tensor.green_tensor_base import GreenTensorBase, get_electric_permitivity
from pairinteraction.units import QuantityScalar, ureg
import green_tensor.utils as utils

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

    def _get_dipole_dipole_au(self, omega_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for a cavity in atomic units.

        Args:
            omega_au: The frequency in atomic units at which to evaluate the Green tensor.

        Returns:
            The dipole dipole Green tensor in cartesian coordinates as a 3x3 array in atomic units.

        """
        if self.pos1_au is None or self.pos2_au is None:
            raise RuntimeError("Atom positions have to be set before calculating the Green tensor.")

        gt = np.zeros((3, 3), dtype=complex)

        au_to_meter = ureg.Quantity(1, "atomic_unit_of_length").to("meter").magnitude
        pos1 = np.array(self.pos1_au) * au_to_meter
        pos2 = np.array(self.pos2_au) * au_to_meter
        epsilon = get_electric_permitivity(self.epsilon, omega_au, "hartree")

        surface1_z = self.surface1_z_au * au_to_meter
        surface1_epsilon = get_electric_permitivity(self.surface1_epsilon, omega_au, "hartree")
        surface2_z = self.surface2_z_au * au_to_meter
        surface2_epsilon = get_electric_permitivity(self.surface2_epsilon, omega_au, "hartree")

        omega = ureg.Quantity(omega_au, "hartree").to("Hz", "spectroscopy").magnitude  # this is the angular frequency

        # TODO calculate Green tensor
#        raise NotImplementedError("GreenTensorFreeSpace is not yet implemented yet.")

        ''' Assumption for the system: The atoms are located at positions pos1 and pos2 at z_A=z_B=h/2. '''

        # If lower surface of the system surface2_z =/ 0, we need to shift the positions accordingly
        pos1_shifted = pos1 - np.array([0, 0, surface2_z])
        pos2_shifted = pos2 - np.array([0, 0, surface2_z])

        # Set height h
        h = surface1_z - surface2_z

        gt = utils.green_tensor_total(pos1_shifted, pos2_shifted, omega, epsilon, surface1_epsilon, surface2_epsilon, h)


        return gt
