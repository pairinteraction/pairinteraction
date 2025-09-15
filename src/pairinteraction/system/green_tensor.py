# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Sequence
from typing import TYPE_CHECKING, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction.units import QuantityArray, QuantityScalar

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.units import Dimension, NDArray, PintArray, PintFloat

    Quantity = TypeVar("Quantity", float, "PintFloat")


class GreenTensor:
    """Green tensor for the multipole pair interactions.

    This class allows to define custom constant or frequency-dependent Green tensors,
    which can then be used for the interaction of a :class:`SystemPair` (see :meth:`SystemPair.set_green_tensor`).

    Examples:
        >>> import pairinteraction as pi
        >>> gt = pi.GreenTensor()
        >>> distance_mum = 5
        >>> tensor = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / distance_mum**3
        >>> tensor_unit = "hartree / (e^2 micrometer^3)"
        >>> gt.set_from_cartesian(1, 1, tensor, tensor_unit)
        GreenTensorReal(...)
        >>> print(gt.get_spherical(1, 1, unit=tensor_unit).diagonal())
        [ 0.008 -0.016  0.008]

    """

    _cpp: _backend.GreenTensorComplex
    _cpp_type = _backend.GreenTensorComplex

    def __init__(self) -> None:
        """Initialize a new Green tensor object.

        The actual tensor can be set afterwards via the :meth:`set_from_cartesian` method.
        """
        self._cpp = self._cpp_type()

    def __repr__(self) -> str:
        return f"{type(self).__name__}(...)"

    def __str__(self) -> str:
        return self.__repr__()

    @overload
    def set_from_cartesian(
        self, kappa1: int, kappa2: int, tensor: "NDArray", tensor_unit: Optional[str] = None
    ) -> "Self": ...

    @overload
    def set_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensor: Sequence["NDArray"],
        tensor_unit: Optional[str],
        omegas: Sequence[float],
        omegas_unit: Optional[str] = None,
    ) -> "Self": ...

    @overload
    def set_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensor: Sequence["NDArray"],
        *,
        omegas: Sequence[float],
        omegas_unit: Optional[str] = None,
    ) -> "Self": ...

    def set_from_cartesian(
        self,
        kappa1: int,
        kappa2: int,
        tensor: Union["NDArray", Sequence["NDArray"]],
        tensor_unit: Optional[str] = None,
        omegas: Optional[Sequence[float]] = None,
        omegas_unit: Optional[str] = None,
    ) -> "Self":
        """Set the entries of the Green tensor.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            tensor: The green tensor in cartesian coordinates.
                Either a single tensor to set a constant green tensor
                or a list of tensors to set a frequency-dependent green tensor.
            tensor_unit: The unit of the tensor.
                Default None, which means that the tensor must be given as pint object.
            omegas: Only needed if a list of tensors is given.
                The frequencies of the tensors.
            omegas_unit: The unit of the frequencies.
                Default None, which means that the frequencies must be given as pint object.


        """
        dimension: list[Dimension] = self._get_unit_dimension(kappa1, kappa2)
        if omegas is None:
            tensor_au = QuantityArray.convert_user_to_au(np.array(tensor), tensor_unit, dimension)
            if tensor_au.shape != (3**kappa1, 3**kappa2) or tensor_au.ndim != 2:
                raise ValueError("The tensor must be a 2D array of shape (3**kappa1, 3**kappa2).")
            self._cpp.create_entries_from_cartesian(kappa1, kappa2, tensor_au)
            return self

        tensors_au = [QuantityArray.convert_user_to_au(np.array(t), tensor_unit, dimension) for t in tensor]
        if not all(t.ndim == 2 for t in tensors_au):
            raise ValueError("The tensor must be a list of 2D arrays.")
        if not all(t.shape == (3**kappa1, 3**kappa2) for t in tensors_au):
            raise ValueError("The tensors must be of shape (3**kappa1, 3**kappa2).")

        omegas_au = QuantityArray.convert_user_to_au(np.array(omegas), omegas_unit, "energy")
        self._cpp.create_entries_from_cartesian(kappa1, kappa2, tensors_au, omegas_au.tolist())
        return self

    @overload
    def get_spherical(
        self,
        kappa1: int,
        kappa2: int,
        omega: Optional[float] = None,
        omega_unit: Optional[str] = None,
        unit: None = None,
    ) -> "PintArray": ...

    @overload
    def get_spherical(
        self, kappa1: int, kappa2: int, omega: Optional[float] = None, omega_unit: Optional[str] = None, *, unit: str
    ) -> "NDArray": ...

    def get_spherical(
        self,
        kappa1: int,
        kappa2: int,
        omega: Optional[float] = None,
        omega_unit: Optional[str] = None,
        unit: Optional[str] = None,
    ) -> Union["PintArray", "NDArray"]:
        """Get the Green tensor in spherical coordinates for the given indices kappa1 and kappa2.

        For kappa == 1 the spherical basis is [p_{1,-1}, p_{1,0}, p_{1,1}].
        For kappa == 2 the spherical basis is [p_{2,-2}, p_{2,-1}, p_{2,0}, p_{2,1}, p_{2,2}, p_{0,0}].

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            omega: The frequency at which to evaluate the Green tensor.
                Only needed if the Green tensor is frequency dependent.
            omega_unit: The unit of the frequency.
                Default None, which means that the frequency must be given as pint object.
            unit: The unit to which to convert the result.
                Default None, which means that the result is returned as pint object.

        Returns:
            The Green tensor as a 2D array.

        """
        entries_cpp = self._cpp.get_spherical_entries(kappa1, kappa2)
        omega_au = QuantityScalar.convert_user_to_au(omega, omega_unit, "energy") if omega is not None else None

        dim1 = 3 if kappa1 == 1 else 6
        dim2 = 3 if kappa2 == 1 else 6
        tensor_au = np.zeros((dim1, dim2))
        for entry_cpp in entries_cpp:
            if isinstance(entry_cpp, (_backend.ConstantEntryReal, _backend.ConstantEntryComplex)):
                val = entry_cpp.val()
            else:
                if omega_au is None:
                    raise ValueError("The Green tensor has omega dependent entries, please specify an omega.")
                val = entry_cpp.val(omega_au)
            tensor_au[entry_cpp.row(), entry_cpp.col()] = val

        return QuantityArray.convert_au_to_user(tensor_au, self._get_unit_dimension(kappa1, kappa2), unit)

    def _get_unit_dimension(self, kappa1: int, kappa2: int) -> list["Dimension"]:
        return ["green_tensor_00", *["inverse_distance" for _ in range(kappa1 + kappa2 + 1)]]  # type: ignore [list-item]


class GreenTensorReal(GreenTensor):
    _cpp: _backend.GreenTensorReal  # type: ignore [assignment]
    _cpp_type = _backend.GreenTensorReal  # type: ignore [assignment]
