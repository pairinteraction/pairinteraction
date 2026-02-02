# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Callable, overload

import numpy as np

from pairinteraction.green_tensor.green_tensor_interpolator import GreenTensorInterpolator
from pairinteraction.units import QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.units import (
        ArrayLike,
        NDArray,
        PintArray,  # needed for sphinx to recognize PintArrayLike
        PintArrayLike,
        PintFloat,
    )


class GreenTensorBase(ABC):
    pos1_au: NDArray | None
    pos2_au: NDArray | None
    epsilon: complex | Callable[[PintFloat], complex]

    def __init__(self) -> None:
        self.pos1_au = None
        self.pos2_au = None
        self.epsilon = 1.0

    def set_atom_positions(
        self, pos1: ArrayLike | PintArrayLike, pos2: ArrayLike | PintArrayLike, unit: str | None = None
    ) -> Self:
        """Set the positions of the two interacting atoms.

        Args:
            pos1: Position of the first atom in the given unit.
            pos2: Position of the second atom in the given unit.
            unit: The unit of the distance, e.g. "micrometer".
                Default None expects a `pint.Quantity`.

        """
        self.pos1_au = np.array([QuantityScalar.convert_user_to_au(v, unit, "distance") for v in pos1])
        self.pos2_au = np.array([QuantityScalar.convert_user_to_au(v, unit, "distance") for v in pos2])
        return self

    def set_electric_permitivity(self, epsilon: complex | Callable[[PintFloat], complex]) -> Self:
        """Set the electric permittivity for the space between the two atoms.

        By default, the electric permittivity is set to 1 (vacuum).

        Args:
            epsilon: The electric permittivity (dimensionless).

        """
        self.epsilon = epsilon
        return self

    @overload
    def get(
        self,
        kappa1: int,
        kappa2: int,
        omega: float = 0,
        omega_unit: str | None = None,
        unit: None = None,
    ) -> PintArray: ...

    @overload
    def get(
        self, kappa1: int, kappa2: int, omega: float = 0, omega_unit: str | None = None, *, unit: str
    ) -> NDArray: ...

    def get(
        self,
        kappa1: int,
        kappa2: int,
        omega: float = 0,
        omega_unit: str | None = None,
        unit: str | None = None,
    ) -> PintArray | NDArray:
        """Calculate the Green tensor in cartesian coordinates for the given indices kappa1 and kappa2.

        kappa == 1 corresponds to dipole operator (-> the basis is [x, y, z]),

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            omega: The angular frequency at which to evaluate the Green tensor.
                Only needed if the Green tensor is frequency dependent.
            omega_unit: The unit of the angular frequency.
                Default None, which means that the angular frequency must be given as pint object.
            unit: The unit to which to convert the result.
                Default None, which means that the result is returned as pint object.

        Returns:
            The Green tensor as a 2D array.

        """
        if kappa1 == 1 and kappa2 == 1:
            return self.get_dipole_dipole(omega, omega_unit, unit=unit)
        raise NotImplementedError("Only dipole-dipole Green tensors are currently implemented.")

    @overload
    def get_dipole_dipole(
        self, omega: float | PintFloat, omega_unit: str | None = None, unit: None = None
    ) -> PintArray: ...

    @overload
    def get_dipole_dipole(self, omega: float | PintFloat, omega_unit: str | None = None, *, unit: str) -> NDArray: ...

    def get_dipole_dipole(
        self, omega: float | PintFloat, omega_unit: str | None = None, unit: str | None = None
    ) -> PintArray | NDArray:
        """Calculate the dipole dipole Green tensor in cartesian coordinates.

        This is a 3x3 matrix in the basis [x1, y1, z1] x [x2, y2, z2].

        Args:
            omega: The angular frequency at which to evaluate the Green tensor.
                Only needed if the Green tensor is frequency dependent.
            omega_unit: The unit of the angular frequency.
                Default None, which means that the angular frequency must be given as pint object.
            unit: The unit to which to convert the green tensor (e.g. "1/m").
                Default None, which means that the result is returned as pint object.

        Returns:
            The dipole dipole Green tensor in cartesian coordinates as a 3x3 array.

        """
        omega_au = QuantityScalar.convert_user_to_au(omega, omega_unit, "energy")
        gt_au = self._get_dipole_dipole_au(omega_au)
        return QuantityArray.convert_au_to_user(gt_au, "green_tensor_dd", unit)

    @abstractmethod
    def _get_dipole_dipole_au(self, omega_au: float) -> NDArray: ...

    def get_green_tensor_interpolator(
        self, omega_min: float, omega_max: float, omega_steps: int, omega_unit: str
    ) -> GreenTensorInterpolator:
        """Get a GreenTensorInterpolator from this Green tensor.

        The GreenTensorInterpolator can be used for the interaction of a SystemPair.

        Returns:
            A GreenTensorInterpolator that interpolates the Green tensor at given frequency points.

        """
        omega_list = list(np.linspace(omega_min, omega_max, omega_steps))
        gt_list = [self.get_dipole_dipole(omega, omega_unit, unit="hartree") for omega in omega_list]

        gti = GreenTensorInterpolator()
        gti.set_list_from_cartesian(1, 1, gt_list, omega_list, tensors_unit="hartree", omegas_unit=omega_unit)

        return gti


def get_electric_permitivity(
    epsilon: complex | Callable[[PintFloat], complex], omega: float, omega_unit: str | None = None
) -> complex:
    """Get the electric permittivity for the given frequency.

    Args:
        epsilon: The electric permittivity (dimensionless) or a callable function that returns it.
        omega: The angular frequency at which to evaluate the permittivity.
            Only needed if the permittivity is frequency dependent.
        omega_unit: The unit of the angular frequency.
            Default None, which means that the angular frequency must be given as pint object.

    Returns:
        The electric permittivity at the given angular frequency.

    """
    omega_au = QuantityScalar.convert_user_to_au(omega, omega_unit, "energy")
    if np.isscalar(epsilon):
        return epsilon  # type: ignore [return-value]
    if callable(epsilon):
        return epsilon(ureg.Quantity(omega_au, "hartree"))
    raise TypeError("epsilon must be either a complex number or a callable function.")
