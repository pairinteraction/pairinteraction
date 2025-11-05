# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase

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
        raise NotImplementedError("GreenTensorFreeSpace is not yet implemented yet.")

        # TODO
