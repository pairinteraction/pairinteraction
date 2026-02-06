# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const

from pairinteraction.green_tensor import utils
from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase, get_electric_permittivity
from pairinteraction.units import QuantityScalar, ureg

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.green_tensor.green_tensor_base import Permitivity
    from pairinteraction.units import (
        ArrayLike,
        NDArray,
        PintArray,  # noqa: F401  # required for sphinx
        PintArrayLike,
        PintFloat,
    )


class GreenTensorSurface(GreenTensorBase):
    """Green tensor for two atoms near a single infinite surface.

    Examples:
        >>> from pairinteraction.green_tensor import GreenTensorSurface
        >>> gt = GreenTensorSurface([0, 0, 0], [10, 0, 0], z=-5, unit="micrometer")
        >>> transition_energy = 2  # h * GHz
        >>> gt_dipole_dipole = gt.get(1, 1, transition_energy, "planck_constant * GHz")
        >>> print(f"{gt_dipole_dipole[0, 0]:.2f}")
        172.51 / bohr

    """

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
        If not specified otherwise (see `set_relative_permittivities`), the surface is treated as a perfect mirror.


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
        self.surface_epsilon: Permitivity = 1e9  # Almost perfect mirror # TODO make utils be able to handle inf

    def set_relative_permittivities(self, epsilon: Permitivity, surface_epsilon: Permitivity) -> Self:
        """Set the relative permittivities of the system.

        Args:
            epsilon: The relative permittivity (dimensionless) of the medium inside the cavity.
            surface_epsilon: The relative permittivity (dimensionless) of the surface.


        """
        self.epsilon = epsilon
        self.surface_epsilon = surface_epsilon
        return self

    def _get_scaled_dipole_dipole_au(self, transition_energy_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for a single surface in atomic units.

        Args:
            transition_energy_au: The transition energy in atomic units at which to evaluate the Green tensor.

        Returns:
            The dipole dipole Green tensor in cartesian coordinates as a 3x3 array in atomic units (i.e. 1/bohr).

        """
        au_to_meter: float = ureg.Quantity(1, "atomic_unit_of_length").to("meter").magnitude
        pos1_m = np.array(self.pos1_au) * au_to_meter
        pos2_m = np.array(self.pos2_au) * au_to_meter
        epsilon = get_electric_permittivity(self.epsilon, transition_energy_au, "hartree")

        omega_hz = ureg.Quantity(transition_energy_au, "hartree").to("hbar Hz", "spectroscopy").magnitude

        # Assume two surfaces, where the further apart atom is located in the center
        # but the second surface has the same permittivity as the inbetween medium
        z1_m = self.surface_z_au * au_to_meter
        height = 2 * max(abs(pos1_m[2] - z1_m), abs(pos2_m[2] - z1_m))
        if pos1_m[2] < z1_m and pos2_m[2] < z1_m:
            z2_m = z1_m - height
        elif pos1_m[2] > z1_m and pos2_m[2] > z1_m:
            z2_m = z1_m + height
        else:
            raise ValueError("Both atoms must be located either above or below the surface.")

        epsilon1 = get_electric_permittivity(self.surface_epsilon, transition_energy_au, "hartree")
        epsilon2 = epsilon

        # unit: # m^(-3) [hbar]^(-1) [epsilon_0]^(-1)
        gt = utils.green_tensor_total(
            pos1_m, pos2_m, z1_m, z2_m, omega_hz, epsilon, epsilon1, epsilon2, only_real_part=True
        )
        to_au = au_to_meter ** (-3) * ((4 * np.pi) ** (-1)) / (const.epsilon_0 * const.hbar)
        # hbar * epsilon_0 = (4*np.pi)**(-1) in atomic units
        return np.real(gt) / to_au
