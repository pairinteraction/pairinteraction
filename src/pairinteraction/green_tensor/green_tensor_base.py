# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Callable, overload

import numpy as np

from pairinteraction.units import QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    from collections.abc import Collection

    from typing_extensions import TypeAlias

    from pairinteraction.green_tensor.green_tensor_interpolator import GreenTensorInterpolator
    from pairinteraction.units import (
        ArrayLike,
        Dimension,
        NDArray,
        PintArray,  # needed for sphinx to recognize PintArrayLike
        PintArrayLike,
        PintFloat,
    )

    PermittivityLike: TypeAlias = complex | Callable[[PintFloat], complex]


class GreenTensorBase(ABC):
    epsilon: complex | Callable[[PintFloat], complex]

    def __init__(
        self,
        pos1: ArrayLike | PintArrayLike,
        pos2: ArrayLike | PintArrayLike,
        unit: str | None = None,
        static_limit: bool = False,
        interaction_order: int = 3,
    ) -> None:
        self.pos1_au = np.array([QuantityScalar.convert_user_to_au(v, unit, "distance") for v in pos1])
        self.pos2_au = np.array([QuantityScalar.convert_user_to_au(v, unit, "distance") for v in pos2])
        self.static_limit = static_limit

        self.interaction_order = interaction_order
        if interaction_order != 3:
            raise NotImplementedError(
                "Only interaction order 3 (dipole-dipole) is currently implemented for custom Green tensors."
            )

        self.epsilon = 1.0

    @overload
    def get(
        self,
        kappa1: int,
        kappa2: int,
        transition_energy: float | PintFloat,
        transition_energy_unit: str | None = None,
        unit: None = None,
        *,
        scaled: bool = False,
    ) -> PintArray: ...

    @overload
    def get(
        self,
        kappa1: int,
        kappa2: int,
        transition_energy: float,
        transition_energy_unit: str,
        unit: str,
        *,
        scaled: bool = False,
    ) -> NDArray: ...

    @overload
    def get(
        self,
        kappa1: int,
        kappa2: int,
        transition_energy: PintFloat,
        *,
        unit: str,
        scaled: bool = False,
    ) -> NDArray: ...

    def get(
        self,
        kappa1: int,
        kappa2: int,
        transition_energy: float | PintFloat,
        transition_energy_unit: str | None = None,
        unit: str | None = None,
        *,
        scaled: bool = False,
    ) -> PintArray | NDArray:
        """Calculate the Green tensor in cartesian coordinates for the given ranks kappa1, kappa2 and frequency omega.

        kappa = 1 corresponds to dipole operator with the cartesian basis: [x, y, z]

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            transition_energy: The transition energy at which to evaluate the Green tensor.
                Only needed if the Green tensor is frequency dependent.
            transition_energy_unit: The unit of the transition energy.
                Default None, which means that the transition energy must be given as pint object.
            unit: The unit to which to convert the result.
                Default None, which means that the result is returned as pint object.
            scaled: If True, the Green tensor is returned with the prefactor for the interaction
                already included (the unit has to be adopted accordingly).
                Default False, which means that the bare Green tensor is returned.

        Returns:
            The Green tensor as a 2D array in cartesian coordinates.

        """
        omega_au = QuantityScalar.convert_user_to_au(transition_energy, transition_energy_unit, "energy")
        omega_for_calculation = 0 if self.static_limit else omega_au

        if kappa1 == 1 and kappa2 == 1:
            scaled_gt_au = self._get_scaled_dipole_dipole_au(omega_for_calculation)
        else:
            raise NotImplementedError("Only dipole-dipole Green tensors are currently implemented.")

        prefactor = 1 if scaled else self._get_prefactor_au(kappa1, kappa2, omega_au)
        dimension = self._get_dimension(kappa1, kappa2, scaled)
        return QuantityArray.convert_au_to_user(scaled_gt_au / prefactor, dimension, unit)

    @abstractmethod
    def _get_scaled_dipole_dipole_au(self, transition_energy_au: float) -> NDArray: ...

    @staticmethod
    def _get_prefactor_au(kappa1: int, kappa2: int, transition_energy_au: float) -> float:
        r"""Get the prefactor to get the interaction strength from the Green tensor.

        The interaction between two dipole moments is given as (see e.g. https://arxiv.org/pdf/2303.13564)
        .. math::
            V_{\alpha\beta} = \frac{\omega^2}{\hbar \epsilon_0 c^2}
                d_\alpha^T \mathrm{Re}\{G(r_\alpha, r_\beta, \omega)\} d_\beta

        This functions returns the prefactor
        .. math::
            \frac{\omega^2}{\hbar \epsilon_0 c^2}

        In C++ we use the convention that the Green tensor already contains this prefactor.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            transition_energy_au: The transition energy at which to evaluate the Green tensor in atomic units (hartree).

        """
        if kappa1 == 1 and kappa2 == 1:
            speed_of_light_au: float = ureg.Quantity(1, "speed_of_light").to_base_units().m
            factor = (transition_energy_au / speed_of_light_au) ** 2
            factor /= 1 / (4 * np.pi)  # hbar * epsilon_0 = 1 / (4*np.pi) in atomic units
            factor /= (2 * np.pi) ** 2  # TODO check omega angular or normal frequency mistake somewhere?
            return factor
        raise NotImplementedError("Only dipole-dipole Green tensor prefactor is currently implemented.")

    @staticmethod
    def _get_dimension(kappa1: int, kappa2: int, scaled: bool) -> list[Dimension]:
        dimension_dd: Dimension = "scaled_green_tensor_dd" if scaled else "green_tensor_dd"
        dimension: list[Dimension] = [dimension_dd] + ["inverse_distance"] * (kappa1 + kappa2 - 2)
        return dimension

    @overload
    def get_interpolator(self, *, use_real: bool) -> GreenTensorInterpolator: ...

    @overload
    def get_interpolator(
        self,
        transition_energies: Collection[PintFloat] | PintArray,
        transition_energies_unit: None = None,
        *,
        use_real: bool,
    ) -> GreenTensorInterpolator: ...

    @overload
    def get_interpolator(
        self, transition_energies: Collection[float] | NDArray, transition_energies_unit: str, *, use_real: bool
    ) -> GreenTensorInterpolator: ...

    def get_interpolator(
        self,
        transition_energies: Collection[PintFloat] | PintArray | Collection[float] | NDArray | None = None,
        transition_energies_unit: str | None = None,
        *,
        use_real: bool,
    ) -> GreenTensorInterpolator:
        """Get a GreenTensorInterpolator from this Green tensor.

        The GreenTensorInterpolator can be used for the interaction of a SystemPair.

        Returns:
            A GreenTensorInterpolator that interpolates the Green tensor at given frequency points.

        """
        GTIClass: type[GreenTensorInterpolator]  # noqa: N806
        if use_real:
            from pairinteraction.green_tensor.green_tensor_interpolator import GreenTensorInterpolatorReal as GTIClass
        else:
            from pairinteraction.green_tensor.green_tensor_interpolator import GreenTensorInterpolator as GTIClass

        if self.static_limit and not (transition_energies is None and transition_energies_unit is None):
            raise ValueError("You must not specify a frequency range when static limit is set.")

        if self.static_limit:
            gti = GTIClass()
            scaled_gt_au = self.get(1, 1, 0, scaled=True)
            gti.set_constant(1, 1, scaled_gt_au)
            return gti

        if transition_energies is not None:
            omegas_pint = [
                QuantityScalar.convert_user_to_pint(omega, transition_energies_unit, "energy")
                for omega in transition_energies
            ]
            gti = GTIClass()
            scaled_gt_list = [self.get(1, 1, omega, scaled=True) for omega in omegas_pint]
            gti.set_list(1, 1, scaled_gt_list, omegas_pint, from_scaled=True)
            return gti

        raise ValueError(
            "You must either specify transition_energies or set static_limit to True to get an interpolator object."
        )


def evaluate_relative_permittivity(
    epsilon: PermittivityLike, transition_energy: float, transition_energy_unit: str | None = None
) -> complex:
    """Get the electric permittivity for the given frequency.

    Args:
        epsilon: The electric permittivity (dimensionless) or a callable function that returns it.
        transition_energy: The angular frequency at which to evaluate the permittivity.
            Only needed if the permittivity is frequency dependent.
        transition_energy_unit: The unit of the angular frequency.
            Default None, which means that the angular frequency must be given as pint object.

    Returns:
        The electric permittivity at the given angular frequency.

    """
    omega_au = QuantityScalar.convert_user_to_au(transition_energy, transition_energy_unit, "energy")
    if np.isscalar(epsilon):
        return epsilon  # type: ignore [return-value]
    if callable(epsilon):
        return epsilon(ureg.Quantity(omega_au, "hartree"))
    raise TypeError("epsilon must be either a complex number or a callable function.")
