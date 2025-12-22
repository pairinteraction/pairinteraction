# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from green_tensor import utils
from green_tensor.green_tensor_base import GreenTensorBase, get_electric_permitivity
from pairinteraction.constants import meter_to_au

from pairinteraction.units import ureg

if TYPE_CHECKING:
    from pairinteraction.units import NDArray


class GreenTensorFreeSpace(GreenTensorBase):
    def _get_dipole_dipole_au(self, omega_au: float) -> NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates for free space in atomic units.

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
        omega = ureg.Quantity(omega_au, "hartree").to("Hz", "spectroscopy").magnitude  # this is the angular frequency

        gt = utils.green_tensor_homogeneous(pos1, pos2, omega, epsilon)
        return gt * 1 / (meter_to_au)
