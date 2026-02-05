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


class GreenTensorCavity(GreenTensorBase):
    """Green tensor for two atoms in a cavity (between two infinite planar surfaces).

    Examples:
        >>> from pairinteraction.green_tensor import GreenTensorCavity
        >>> gt = GreenTensorCavity([0, 0, 0], [10, 0, 0], z1=-5, z2=5, unit="micrometer")
        >>> transition_energy = 2  # h * GHz
        >>> gt_dipole_dipole = gt.get(1, 1, transition_energy, "planck_constant * GHz")
        >>> print(f"{gt_dipole_dipole[0, 0]:.2f}")
        151.77 / bohr

    """

    def __init__(
        self,
        pos1: ArrayLike | PintArrayLike,
        pos2: ArrayLike | PintArrayLike,
        z1: float | PintFloat,
        z2: float | PintFloat,
        unit: str | None = None,
        static_limit: bool = False,
        interaction_order: int = 3,
    ) -> None:
        """Create a Green tensor for two atoms inside a planar cavity formed by two infinite surfaces.

        The two surfaces of the cavity are assumed to be infinite in the x-y plane.
        If not specified otherwise (see `set_relative_permittivities`), the surfaces are treated as perfect mirrors.

        Args:
            pos1: Position of the first atom in the given unit.
            pos2: Position of the second atom in the given unit.
            z1: The z-position of the first surface in the given unit.
            z2: The z-position of the second surface in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.
            static_limit: If True, the static limit is used.
                Default False.
            interaction_order: The order of interaction, e.g., 3 for dipole-dipole.
                Defaults to 3.

        """
        super().__init__(pos1, pos2, unit, static_limit, interaction_order)
        self.surface1_z_au = QuantityScalar.convert_user_to_au(z1, unit, "distance")
        self.surface2_z_au = QuantityScalar.convert_user_to_au(z2, unit, "distance")
        # Almost perfect mirrors # TODO make utils be able to handle inf
        self.surface1_epsilon: Permitivity = 1e9
        self.surface2_epsilon: Permitivity = 1e9

    def set_relative_permittivities(self, epsilon: Permitivity, epsilon1: Permitivity, epsilon2: Permitivity) -> Self:
        """Set the relative permittivities of the system.

        Args:
            epsilon: The relative permittivity (dimensionless) of the medium inside the cavity.
            epsilon1: The relative permittivity (dimensionless) of the first surface.
            epsilon2: The relative permittivity (dimensionless) of the second surface.


        """
        self.epsilon = epsilon
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
        epsilon = get_electric_permittivity(self.epsilon, omega_au, "hartree")

        omega = ureg.Quantity(omega_au, "hartree").to("hbar Hz")
        omega_hz = omega.magnitude

        surface1_z_m = self.surface1_z_au * au_to_meter
        surface2_z_m = self.surface2_z_au * au_to_meter

        height = abs(surface1_z_m - surface2_z_m)
        pos1_shifted_m = pos1_m - np.array([0, 0, min(surface1_z_m, surface2_z_m)])
        pos2_shifted_m = pos2_m - np.array([0, 0, min(surface1_z_m, surface2_z_m)])

        if not 0 < pos1_shifted_m[2] < height or not 0 < pos2_shifted_m[2] < height:
            raise ValueError("Both atoms must be located inside the cavity between the two surfaces.")

        if surface2_z_m < surface1_z_m:
            epsilon_top = get_electric_permittivity(self.surface1_epsilon, omega_au, "hartree")
            epsilon_bottom = get_electric_permittivity(self.surface2_epsilon, omega_au, "hartree")
        else:
            epsilon_top = get_electric_permittivity(self.surface2_epsilon, omega_au, "hartree")
            epsilon_bottom = get_electric_permittivity(self.surface1_epsilon, omega_au, "hartree")

        # unit: # m^(-3) [hbar]^(-1) [epsilon_0]^(-1)
        gt = utils.green_tensor_total(
            pos1_shifted_m, pos2_shifted_m, omega_hz, epsilon, epsilon_top, epsilon_bottom, height, only_real_part=True
        )
        to_au = au_to_meter ** (-3) * ((4 * np.pi) ** (-1)) / (const.epsilon_0 * const.hbar)
        # hbar * epsilon_0 = (4*np.pi)**(-1) in atomc units
        return np.real(gt) / to_au
