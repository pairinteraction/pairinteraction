# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from pairinteraction.green_tensor import utils
from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase, get_electric_permitivity
from pairinteraction.units import QuantityScalar, ureg

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

        au_to_meter = ureg.Quantity(1, "atomic_unit_of_length").to("meter").magnitude
        pos1_m = np.array(self.pos1_au) * au_to_meter
        pos2_m = np.array(self.pos2_au) * au_to_meter
        epsilon = get_electric_permitivity(self.epsilon, omega_au, "hartree")

        omega = ureg.Quantity(omega_au, "hartree").to("Hz", "spectroscopy")  # this is the angular frequency
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

        gt = utils.green_tensor_total(
            pos1_shifted_m, pos2_shifted_m, omega_hz, epsilon, epsilon_top, epsilon_bottom, height, only_real_part=True
        )

        # see green_tensor_free_space for details on the prefactors
        gt *= 4 * np.pi * au_to_meter * (omega / ureg.speed_of_light).to_base_units().m ** 2
        return np.real(gt)
