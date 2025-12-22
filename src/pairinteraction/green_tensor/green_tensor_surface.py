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


class GreenTensorSurface(GreenTensorBase):
    def __init__(self) -> None:
        super().__init__()
        self.surface_z_au = 0.0

    def set_surface_position(self, z: float, unit: str | None = None) -> Self:
        """Set the position of the surface along the z-axis.

        The surface is assumed to be infinite in the x-y plane.
        By default (i.e. if this method was not called), the surface is located at z=0.

        Args:
            z: The z-position of the surface in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        self.surface_z_au = QuantityScalar.convert_user_to_au(z, unit, "distance")
        return self

    def set_electric_permitivity_surface(self, epsilon: complex | Callable[[PintFloat], complex]) -> Self:
        """Set the electric permittivity for the surface.

        Args:
            epsilon: The electric permittivity (dimensionless).

        """
        self.surface_epsilon = epsilon
        return self


    def _get_dipole_dipole_au(self, omega_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for a single surface in atomic units.

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

        surface_z = self.surface_z_au * au_to_meter
        surface_epsilon = get_electric_permitivity(self.surface_epsilon, omega_au, "hartree")

        omega = ureg.Quantity(omega_au, "hartree").to("Hz", "spectroscopy").magnitude  # this is the angular frequency

        # TODO calculate Green tensor
#        raise NotImplementedError("GreenTensorFreeSpace is not yet implemented yet.")

        ''' Assumption for the system: The atoms are located at positions pos1 and pos2 above a single surface at z=surface_z.
            It should be pos1[2], pos2[2] > surface_z for this to be valid.'''

        # If z_surface =/ 0, we need to shift the positions accordingly
        pos1_shifted = pos1 - np.array([0, 0, surface_z])
        pos2_shifted = pos2 - np.array([0, 0, surface_z])

        # Set height h to 2 * z_A
        h = 2 * pos1_shifted[2]

        # Set the permittivities of the system
        epsilon0 = epsilon  # permittivity of the medium above the surface
        epsilon1 = 1.0 # upper layer is vacuum, since we only have one surface
        epsilon2 = surface_epsilon  # permittivity of the surface

        gt = utils.green_tensor_total(pos1_shifted, pos2_shifted, omega, epsilon0, epsilon1, epsilon2, h)

        return gt