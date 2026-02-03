# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const

from pairinteraction.green_tensor import utils
from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase, get_electric_permitivity
from pairinteraction.units import ureg

if TYPE_CHECKING:
    from pairinteraction.units import NDArray


class GreenTensorFreeSpace(GreenTensorBase):
    def _get_scaled_dipole_dipole_au(self, omega_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for free space in atomic units.

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

        # unit: # m^(-3) [hbar]^(-1) [epsilon_0]^(-1)
        gt = utils.green_tensor_homogeneous(pos1_m, pos2_m, omega_hz, epsilon, only_real_part=True)
        to_au = ureg.Quantity(1, "m**3") * ((4 * np.pi) ** (-1)) / (const.epsilon_0 * const.hbar)
        # hbar * epsilon_0 = (4*np.pi)**(-1) in atomc units
        return np.real(gt) / to_au.to_base_units().magnitude
