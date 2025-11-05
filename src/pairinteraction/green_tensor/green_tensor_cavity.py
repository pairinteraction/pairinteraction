# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase
from pairinteraction.units import QuantityScalar

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
        self, epsilon1: float | Callable[[PintFloat], float], epsilon2: float | Callable[[PintFloat], float]
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
        raise NotImplementedError("GreenTensorCavity is not yet implemented yet.")

        # TODO
