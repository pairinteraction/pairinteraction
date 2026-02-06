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

    from pairinteraction.units import NDArray, PintArray, PintFloat

Coordinates = Literal["cartesian", "spherical"]


class GreenTensorInterpolator:
    """Green tensor interpolator for the multipole pair interactions.

    This class allows to define constant or frequency-dependent Green tensor interpolators,
    which can then be used for the interaction of a :class:`SystemPair`
    (see :meth:`SystemPair._set_green_tensor_interpolator`).

    Examples:
        >>> from pairinteraction.green_tensor import GreenTensorInterpolator
        >>> gt = GreenTensorInterpolator()
        >>> distance_mum = 5
        >>> transition_energy = 1  # planck_constant * GHz
        >>> tensor = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / (distance_mum**3)
        >>> tensor_unit = "hartree / (e * micrometer)^2"
        >>> gt.set_constant(1, 1, tensor, tensor_unit)
        GreenTensorInterpolator(...)
        >>> tensor_sph = gt.get(1, 1, transition_energy, "planck_constant * GHz", unit=tensor_unit, scaled=True)
        >>> print(tensor_sph.diagonal())
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
                (see :meth:`GreenTensorBase._get_prefactor_au`).
            tensor_unit: The unit of the tensor.
                Default None, which means that the tensor must be given as pint object.
            from_coordinates: The coordinate system in which the tensor is given.
                Default "cartesian".

        """
        if from_coordinates != "cartesian":
            raise NotImplementedError("Only cartesian coordinates are currently implemented for set_constant.")
        if tensor.shape != (3**kappa1, 3**kappa2) or tensor.ndim != 2:  # type: ignore [union-attr]
            raise ValueError("The tensor must be a 2D array of shape (3**kappa1, 3**kappa2).")

        dimension = GreenTensorBase._get_dimension(kappa1, kappa2, scaled=True)
        scaled_tensor_au = QuantityArray.convert_user_to_au(tensor, tensor_unit, dimension)
        self._cpp.create_entries_from_cartesian(kappa1, kappa2, scaled_tensor_au)
        return self

    @overload
    def set_list(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[PintArray],
        transition_energies: Collection[PintFloat],
        tensors_unit: None = None,
        transition_energies_unit: None = None,
        *,
        from_coordinates: Coordinates = "cartesian",
        from_scaled: bool = False,
    ) -> Self: ...

    @overload
    def set_list(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[NDArray],
        transition_energies: Collection[float],
        tensors_unit: str,
        transition_energies_unit: str,
        *,
        from_coordinates: Coordinates = "cartesian",
        from_scaled: bool = False,
    ) -> Self: ...

    def set_list(
        self,
        kappa1: int,
        kappa2: int,
        tensors: Collection[PintArray] | Collection[NDArray],
        transition_energies: Collection[PintFloat] | Collection[float],
        tensors_unit: str | None = None,
        transition_energies_unit: str | None = None,
        *,
        from_coordinates: Coordinates = "cartesian",
        from_scaled: bool = False,
    ) -> Self:
        """Set the entries of the Green tensor for specified transition energies.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            tensors: A list of frequency-dependent green tensors in cartesian coordinates.
            transition_energies: A list of transition energies at which the green tensors are defined.
            tensors_unit: The unit of the tensor.
                Default None, which means that the tensor must be given as pint object.
            transition_energies_unit: The unit of the transition energies.
                Default None, which means that the transition energies must be given as pint object.
            from_coordinates: The coordinate system in which the tensor is given.
                Default "cartesian".
            from_scaled: Whether the prefactor for the interaction strength
                (see :meth:`GreenTensorBase._get_prefactor_au`) is already included in the given tensor.
                The unit has to be adjusted accordingly. Default False.

        """
        if from_coordinates != "cartesian":
            raise NotImplementedError("Only cartesian coordinates are currently implemented for set_list.")
        if not all(t.ndim == 2 for t in tensors):
            raise ValueError("The tensor must be a list of 2D arrays.")
        if not all(t.shape == (3**kappa1, 3**kappa2) for t in tensors):  # type: ignore [union-attr]
            raise ValueError("The tensors must be of shape (3**kappa1, 3**kappa2).")

        omegas_au = [
            QuantityScalar.convert_user_to_au(omega, transition_energies_unit, "energy")
            for omega in transition_energies
        ]
        prefactors = [1.0] * len(tensors)
        if not from_scaled:
            prefactors = [GreenTensorBase._get_prefactor_au(kappa1, kappa2, omega) for omega in omegas_au]
        dimension = GreenTensorBase._get_dimension(kappa1, kappa2, scaled=from_scaled)
        tensors_au = [QuantityArray.convert_user_to_au(t, tensors_unit, dimension) for t in tensors]
        scaled_tensors_au = [prefactor * t for prefactor, t in zip(prefactors, tensors_au)]
        self._cpp.create_entries_from_cartesian(kappa1, kappa2, scaled_tensors_au, omegas_au)  # type: ignore [arg-type]
        return self

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
        coordinates: Coordinates = "spherical",
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
        coordinates: Coordinates = "spherical",
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
        coordinates: Coordinates = "spherical",
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
        coordinates: Coordinates = "spherical",
    ) -> PintArray | NDArray:
        """Get the Green tensor in the given coordinates for the given ranks kappa1, kappa2 and transition energy.

        kappa = 1 corresponds to dipole operator with the basis
            - spherical: [p_{1,-1}, p_{1,0}, p_{1,1}]
        kappa = 2 corresponds to quadrupole operator with the basis
            - spherical: [p_{2,-2}, p_{2,-1}, p_{2,0}, p_{2,1}, p_{2,2}, p_{0,0}]

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            transition_energy: The transition energy at which to evaluate the Green tensor.
                Use transition_energy=0 for the static limit.
            transition_energy_unit: The unit of the transition energy.
                Default None, which means that the transition energy must be given as pint object (or is 0).
            unit: The unit to which to convert the result.
                Default None, which means that the result is returned as pint object.
            scaled: If True, the Green tensor is returned with the prefactor for the interaction
                already included (the unit has to be adopted accordingly).
                Default False returns the bare Green tensor.
            coordinates: The coordinate system in which to return the tensor.
                Default "spherical".

        Returns:
            The Green tensor as a 2D array.

        """
        if coordinates != "spherical":
            raise NotImplementedError("Only spherical coordinates are currently implemented for get.")

        omega_au = QuantityScalar.convert_user_to_au(transition_energy, transition_energy_unit, "energy")

        entries_cpp = self._cpp.get_spherical_entries(kappa1, kappa2)
        kappa_to_dim = {1: 3, 2: 6}
        dim1, dim2 = kappa_to_dim[kappa1], kappa_to_dim[kappa2]
        tensor_au = np.zeros((dim1, dim2), dtype=complex)
        for entry_cpp in entries_cpp:
            if isinstance(entry_cpp, (_backend.ConstantEntryReal, _backend.ConstantEntryComplex)):
                val = entry_cpp.val()
            else:
                val = entry_cpp.val(omega_au)
            tensor_au[entry_cpp.row(), entry_cpp.col()] = val
        tensor_au = np.real_if_close(tensor_au)

        prefactor = 1 if scaled else GreenTensorBase._get_prefactor_au(kappa1, kappa2, omega_au)
        dimension = GreenTensorBase._get_dimension(kappa1, kappa2, scaled=scaled)
        return QuantityArray.convert_au_to_user(tensor_au / prefactor, dimension, unit)


class GreenTensorInterpolatorReal(GreenTensorInterpolator):
    _cpp: _backend.GreenTensorInterpolatorReal  # type: ignore [assignment]
    _cpp_type = _backend.GreenTensorInterpolatorReal  # type: ignore [assignment]
