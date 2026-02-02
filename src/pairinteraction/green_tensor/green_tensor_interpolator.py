# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction.units import QuantityArray, QuantityScalar, ureg

if TYPE_CHECKING:
    from collections.abc import Collection

    from typing_extensions import Self

    from pairinteraction.units import Dimension, NDArray, PintArray, PintFloat


class GreenTensorInterpolator:
    """Green tensor interpolator for the multipole pair interactions.

    This class allows to define constant or frequency-dependent Green tensor interpolators,
    which can then be used for the interaction of a :class:`SystemPair`
    (see :meth:`SystemPair.set_green_tensor_interpolator`).

    Examples:
        >>> import pairinteraction as pi
        >>> gt = pi.GreenTensorInterpolator()
        >>> distance_mum = 5
        >>> omega = 1
        >>> tensor = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / (omega**2 * distance_mum**3)
        >>> tensor_unit = "hartree / (e^2 micrometer^3)"
        >>> gt.set_constant_from_cartesian(1, 1, tensor, omega, tensor_unit, omega_unit="Hz")
        GreenTensorInterpolator(...)
        >>> print(gt.get_spherical(1, 1, unit=tensor_unit).diagonal())
        [ 0.008 -0.016  0.008]

    """

    _cpp: _backend.GreenTensorInterpolatorComplex
    _cpp_type = _backend.GreenTensorInterpolatorComplex

    def __init__(self) -> None:
        """Initialize a new Green tensor interpolator object.

        The actual tensor can be set afterwards via the
        :meth:`set_constant_from_cartesian` or :meth:`set_list_from_cartesian` method.
        """
        self._cpp = self._cpp_type()

    def __repr__(self) -> str:
        return f"{type(self).__name__}(...)"

    def __str__(self) -> str:
        return self.__repr__()

    @overload
    def set_constant_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensor: PintArray,
        omega: PintFloat,
        tensor_unit: None = None,
        omega_unit: None = None,
    ) -> Self: ...

    @overload
    def set_constant_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensor: NDArray,
        omega: float,
        tensor_unit: str,
        omega_unit: str,
    ) -> Self: ...

    def set_constant_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensor: NDArray | PintArray,
        omega: float | PintFloat,
        tensor_unit: str | None = None,
        omega_unit: str | None = None,
    ) -> Self:
        r"""Set the constant entry of the Green tensor for specified omega.

        Constant means, that :math:`\omega^2 G(\omega)` (which is the quantity that enters the interaction)
        is constant and independent of omega.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            tensor: The green tensor in cartesian coordinates.
            omega: The angular frequency at which the green tensor is defined.
            tensor_unit: The unit of the tensor.
                Default None, which means that the tensor must be given as pint object.
            omega_unit: The unit of the angular frequency.
                Default None, which means that the angular frequency must be given as pint object.


        """
        dimension: list[Dimension] = self._get_unit_dimension(kappa1, kappa2)

        tensor_au = QuantityArray.convert_user_to_au(tensor, tensor_unit, dimension)
        if tensor_au.shape != (3**kappa1, 3**kappa2) or tensor_au.ndim != 2:
            raise ValueError("The tensor must be a 2D array of shape (3**kappa1, 3**kappa2).")
        omega_au = QuantityScalar.convert_user_to_au(omega, omega_unit, "energy")

        tensor_cpp: NDArray = self._get_green_tensor_prefactor(omega_au) * tensor_au
        self._cpp.create_entries_from_cartesian(kappa1, kappa2, tensor_cpp)
        return self

    @overload
    def set_list_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[PintArray] | PintArray,
        omegas: Collection[PintFloat] | PintArray,
        tensors_unit: None = None,
        omegas_unit: None = None,
    ) -> Self: ...

    @overload
    def set_list_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[NDArray],
        omegas: Collection[float],
        tensors_unit: str,
        omegas_unit: str,
    ) -> Self: ...

    def set_list_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[PintArray] | PintArray | Collection[NDArray],
        omegas: Collection[PintFloat] | PintArray | Collection[float],
        tensors_unit: str | None = None,
        omegas_unit: str | None = None,
    ) -> Self:
        """Set the entries of the Green tensor for specified omega.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            tensors: A list of frequency-dependent green tensors in cartesian coordinates.
            omegas: A list of angular frequencies at which the green tensors are defined.
            tensors_unit: The unit of the tensor.
                Default None, which means that the tensor must be given as pint object.
            omegas_unit: The unit of the angular frequencies.
                Default None, which means that the angular frequencies must be given as pint object.

        """
        dimension: list[Dimension] = self._get_unit_dimension(kappa1, kappa2)

        tensors_au = np.array([QuantityArray.convert_user_to_au(t, tensors_unit, dimension) for t in tensors])
        if not all(t.ndim == 2 for t in tensors_au):
            raise ValueError("The tensor must be a list of 2D arrays.")
        if not all(t.shape == (3**kappa1, 3**kappa2) for t in tensors_au):
            raise ValueError("The tensors must be of shape (3**kappa1, 3**kappa2).")
        omegas_au = [QuantityScalar.convert_user_to_au(omega, omegas_unit, "energy") for omega in omegas]
        prefactors = [self._get_green_tensor_prefactor(omega) for omega in omegas_au]

        tensors_cpp = [prefactor * tensor for prefactor, tensor in zip(prefactors, tensors_au)]
        self._cpp.create_entries_from_cartesian(kappa1, kappa2, tensors_cpp, omegas_au)
        return self

    @overload
    def get_spherical(
        self,
        kappa1: int,
        kappa2: int,
        omega: float | PintFloat,
        omega_unit: str | None = None,
        unit: None = None,
    ) -> PintArray: ...

    @overload
    def get_spherical(
        self, kappa1: int, kappa2: int, omega: float | PintFloat, omega_unit: str | None = None, *, unit: str
    ) -> NDArray: ...

    def get_spherical(
        self,
        kappa1: int,
        kappa2: int,
        omega: float | PintFloat,
        omega_unit: str | None = None,
        unit: str | None = None,
    ) -> PintArray | NDArray:
        """Get the Green tensor in spherical coordinates for the given indices kappa1, kappa2 and frequency omega.

        For kappa == 1 the spherical basis is [p_{1,-1}, p_{1,0}, p_{1,1}].
        For kappa == 2 the spherical basis is [p_{2,-2}, p_{2,-1}, p_{2,0}, p_{2,1}, p_{2,2}, p_{0,0}].

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            omega: The angular frequency at which to evaluate the Green tensor.
            omega_unit: The unit of the angular frequency.
                Default None, which means that the angular frequency must be given as pint object.
            unit: The unit to which to convert the result.
                Default None, which means that the result is returned as pint object.

        Returns:
            The Green tensor as a 2D array.

        """
        entries_cpp = self._cpp.get_spherical_entries(kappa1, kappa2)
        omega_au = QuantityScalar.convert_user_to_au(omega, omega_unit, "energy")

        dim1 = 3 if kappa1 == 1 else 6
        dim2 = 3 if kappa2 == 1 else 6
        tensor_au = np.zeros((dim1, dim2), dtype=complex)
        for entry_cpp in entries_cpp:
            if isinstance(entry_cpp, (_backend.ConstantEntryReal, _backend.ConstantEntryComplex)):
                val = entry_cpp.val()
            else:
                val = entry_cpp.val(omega_au)
            tensor_au[entry_cpp.row(), entry_cpp.col()] = val
        tensor_au = np.real_if_close(tensor_au)
        tensor_au /= self._get_green_tensor_prefactor(omega_au)

        return QuantityArray.convert_au_to_user(tensor_au, self._get_unit_dimension(kappa1, kappa2), unit)

    def _get_unit_dimension(self, kappa1: int, kappa2: int) -> list[Dimension]:
        return ["green_tensor_00", *["inverse_distance" for _ in range(kappa1 + kappa2)]]  # type: ignore [list-item]

    def _get_green_tensor_prefactor(self, omega_au: float) -> float:
        r"""Get the prefactor to get the interaction strength from the Green tensor.

        The interaction between two dipole moments is given as (see e.g. https://arxiv.org/pdf/2303.13564)
        .. math::
            V_{\alpha\beta} = \frac{\omega^2}{\hbar \epsilon_0 c^2}
                d_\alpha^T \mathrm{Re}\{G(r_\alpha, r_\beta, \omega)\} d_\beta

        This functions returns the prefactor
        .. math::
            \frac{\omega^2}{\hbar \epsilon_0 c^2}

        In C++ we use the convention that the Green tensor already contains this prefactor.
        """
        speed_of_light_au = ureg.Quantity(1, "speed_of_light").to_base_units().m
        factor = (omega_au / speed_of_light_au) ** 2
        factor /= 1 / (4 * np.pi)  # hbar * epsilon_0 = 1 / (4*np.pi) in atomc units
        factor /= (2 * np.pi) ** 2  # TODO check omega angular or normal frequency mistake somewhere?
        return factor


class GreenTensorInterpolatorReal(GreenTensorInterpolator):
    _cpp: _backend.GreenTensorInterpolatorReal  # type: ignore [assignment]
    _cpp_type = _backend.GreenTensorInterpolatorReal  # type: ignore [assignment]
