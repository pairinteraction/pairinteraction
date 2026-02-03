# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

from typing import TYPE_CHECKING, Literal, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction.green_tensor.green_tensor_base import GreenTensorBase
from pairinteraction.units import QuantityArray, QuantityScalar

if TYPE_CHECKING:
    from collections.abc import Collection

    from typing_extensions import Self

    from pairinteraction.units import Dimension, NDArray, PintArray, PintFloat

Coordinates = Literal["cartesian", "spherical"]


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
        >>> tensor_unit = "1 / micrometer"
        >>> gt.set_constant(1, 1, tensor, omega, tensor_unit, omega_unit="Hz")
        GreenTensorInterpolator(...)
        >>> print(gt.get(1, 1, omega, omega_unit="Hz", unit=tensor_unit).diagonal())
        [ 0.008 -0.016  0.008]

    """

    _cpp: _backend.GreenTensorInterpolatorComplex
    _cpp_type = _backend.GreenTensorInterpolatorComplex

    def __init__(self) -> None:
        """Initialize a new Green tensor interpolator object.

        The actual tensor can be set afterwards via the
        :meth:`set_constant` or :meth:`set_list` method.
        """
        self._cpp = self._cpp_type()

    def __repr__(self) -> str:
        return f"{type(self).__name__}(...)"

    def __str__(self) -> str:
        return self.__repr__()

    @overload
    def set_constant(
        self,
        kappa1: int,
        kappa2: int,
        tensor: PintArray,
        tensor_unit: None = None,
        *,
        from_coordinates: Coordinates = "cartesian",
    ) -> Self: ...

    @overload
    def set_constant(
        self,
        kappa1: int,
        kappa2: int,
        tensor: NDArray,
        tensor_unit: str,
        *,
        from_coordinates: Coordinates = "cartesian",
    ) -> Self: ...

    def set_constant(
        self,
        kappa1: int,
        kappa2: int,
        tensor: NDArray | PintArray,
        tensor_unit: str | None = None,
        *,
        from_coordinates: Coordinates = "cartesian",
    ) -> Self:
        r"""Set the scaled Green tensor to a constant entry.

        Constant means, that :math:`\omega^2 G(\omega)` (which is the quantity that enters the interaction)
        is constant and independent of omega.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            tensor: The scaled green tensor including the prefactor for the interaction strength
                (see :meth:`GreenTensorBase._get_green_tensor_prefactor_au`).
            tensor_unit: The unit of the tensor.
                Default None, which means that the tensor must be given as pint object.
            from_coordinates: The coordinate system in which the tensor is given.
                Default "cartesian".

        """
        if from_coordinates != "cartesian":
            raise NotImplementedError("Only cartesian coordinates are currently implemented for set_constant.")

        dimension: list[Dimension] = ["scaled_green_tensor_dd"]
        dimension += ["inverse_distance" for _ in range(kappa1 + kappa2 - 2)]

        scaled_tensor_au = QuantityArray.convert_user_to_au(tensor, tensor_unit, dimension)
        if scaled_tensor_au.shape != (3**kappa1, 3**kappa2) or scaled_tensor_au.ndim != 2:
            raise ValueError("The tensor must be a 2D array of shape (3**kappa1, 3**kappa2).")

        self._cpp.create_entries_from_cartesian(kappa1, kappa2, scaled_tensor_au)
        return self

    @overload
    def set_list(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[PintArray] | PintArray,
        omegas: Collection[PintFloat] | PintArray,
        tensors_unit: None = None,
        omegas_unit: None = None,
        *,
        from_coordinates: Coordinates = "cartesian",
        prefactor_already_included: bool = False,
    ) -> Self: ...

    @overload
    def set_list(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[NDArray],
        omegas: Collection[float],
        tensors_unit: str,
        omegas_unit: str,
        *,
        from_coordinates: Coordinates = "cartesian",
        prefactor_already_included: bool = False,
    ) -> Self: ...

    def set_list(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[PintArray] | PintArray | Collection[NDArray],
        omegas: Collection[PintFloat] | PintArray | Collection[float],
        tensors_unit: str | None = None,
        omegas_unit: str | None = None,
        *,
        from_coordinates: Coordinates = "cartesian",
        prefactor_already_included: bool = False,
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
            from_coordinates: The coordinate system in which the tensor is given.
                Default "cartesian".
            prefactor_already_included: Whether the prefactor for the interaction strength
                (see :meth:`GreenTensorBase._get_green_tensor_prefactor_au`) is already included in the given tensor.
                If True, the unit has to be adjusted accordingly.
                Default False.

        """
        if from_coordinates != "cartesian":
            raise NotImplementedError("Only cartesian coordinates are currently implemented for set_list.")

        dimension: list[Dimension] = ["scaled_green_tensor_dd" if prefactor_already_included else "green_tensor_dd"]
        dimension += ["inverse_distance" for _ in range(kappa1 + kappa2 - 2)]

        tensors_au = [QuantityArray.convert_user_to_au(t, tensors_unit, dimension) for t in tensors]
        if not all(t.ndim == 2 for t in tensors_au):
            raise ValueError("The tensor must be a list of 2D arrays.")
        if not all(t.shape == (3**kappa1, 3**kappa2) for t in tensors_au):
            raise ValueError("The tensors must be of shape (3**kappa1, 3**kappa2).")

        omegas_au = [QuantityScalar.convert_user_to_au(omega, omegas_unit, "energy") for omega in omegas]
        if not prefactor_already_included:
            prefactors = [GreenTensorBase._get_green_tensor_prefactor_au(omega) for omega in omegas_au]
            tensors_au = [prefactor * tensor for prefactor, tensor in zip(prefactors, tensors_au)]

        self._cpp.create_entries_from_cartesian(kappa1, kappa2, tensors_au, omegas_au)
        return self

    @overload
    def get(
        self,
        kappa1: int,
        kappa2: int,
        omega: float | PintFloat,
        omega_unit: str | None = None,
        unit: None = None,
        *,
        scaled_green_tensor: bool = False,
        coordinates: Coordinates = "spherical",
    ) -> PintArray: ...

    @overload
    def get(
        self,
        kappa1: int,
        kappa2: int,
        omega: float | PintFloat,
        omega_unit: str | None = None,
        *,
        unit: str,
        scaled_green_tensor: bool = False,
        coordinates: Coordinates = "spherical",
    ) -> NDArray: ...

    def get(
        self,
        kappa1: int,
        kappa2: int,
        omega: float | PintFloat,
        omega_unit: str | None = None,
        unit: str | None = None,
        *,
        scaled_green_tensor: bool = False,
        coordinates: Coordinates = "spherical",
    ) -> PintArray | NDArray:
        """Get the Green tensor in the given coordinates for the given ranks kappa1, kappa2 and frequency omega.

        kappa = 1 corresponds to dipole operator with the basis
            - spherical: [p_{1,-1}, p_{1,0}, p_{1,1}]
        kappa = 2 corresponds to quadrupole operator with the basis
            - spherical: [p_{2,-2}, p_{2,-1}, p_{2,0}, p_{2,1}, p_{2,2}, p_{0,0}]

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            omega: The angular frequency at which to evaluate the Green tensor.
                Use omega=0 for the static limit.
            omega_unit: The unit of the angular frequency.
                Default None, which means that the angular frequency must be given as pint object (or is 0).
            unit: The unit to which to convert the result.
                Default None, which means that the result is returned as pint object.
            scaled_green_tensor: If True, the Green tensor is returned with the prefactor for the interaction
                already included (the unit has to be adopted accordingly).
                Default False returns the bare Green tensor.
            coordinates: The coordinate system in which to return the tensor.
                Default "spherical".

        Returns:
            The Green tensor as a 2D array.

        """
        if coordinates != "spherical":
            raise NotImplementedError("Only spherical coordinates are currently implemented for get.")

        entries_cpp = self._cpp.get_spherical_entries(kappa1, kappa2)
        omega_au = QuantityScalar.convert_user_to_au(omega, omega_unit, "energy")

        kappa_to_dim = {1: 3, 2: 6}
        dim1 = kappa_to_dim[kappa1]
        dim2 = kappa_to_dim[kappa2]
        tensor_au = np.zeros((dim1, dim2), dtype=complex)
        for entry_cpp in entries_cpp:
            if isinstance(entry_cpp, (_backend.ConstantEntryReal, _backend.ConstantEntryComplex)):
                val = entry_cpp.val()
            else:
                val = entry_cpp.val(omega_au)
            tensor_au[entry_cpp.row(), entry_cpp.col()] = val
        tensor_au = np.real_if_close(tensor_au)

        dimension: list[Dimension] = ["scaled_green_tensor_dd" if scaled_green_tensor else "green_tensor_dd"]
        dimension += ["inverse_distance" for _ in range(kappa1 + kappa2 - 2)]

        if not scaled_green_tensor:
            tensor_au /= GreenTensorBase._get_green_tensor_prefactor_au(omega_au)

        return QuantityArray.convert_au_to_user(tensor_au, dimension, unit)


class GreenTensorInterpolatorReal(GreenTensorInterpolator):
    _cpp: _backend.GreenTensorInterpolatorReal  # type: ignore [assignment]
    _cpp_type = _backend.GreenTensorInterpolatorReal  # type: ignore [assignment]
