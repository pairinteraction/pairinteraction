# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Sequence
from typing import TYPE_CHECKING, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction.units import QuantityArray, QuantityScalar

if TYPE_CHECKING:
    from pairinteraction.units import Dimension, NDArray, PintArray, PintFloat

    Quantity = TypeVar("Quantity", float, "PintFloat")

UnionCPPGreenTensor = Union[_backend.GreenTensorReal, _backend.GreenTensorComplex]
UnionTypeCPPGreenTensor = Union[type[_backend.GreenTensorReal], type[_backend.GreenTensorComplex]]

UnionCPPEntry = Union[
    _backend.ConstantEntryReal,
    _backend.ConstantEntryComplex,
    _backend.OmegaDependentEntryReal,
    _backend.OmegaDependentEntryComplex,
]
UnionCPPConstantEntry = Union[_backend.ConstantEntryReal, _backend.ConstantEntryComplex]
UnionCPPOmegaDependentEntry = Union[_backend.OmegaDependentEntryReal, _backend.OmegaDependentEntryComplex]


class GreenTensor:
    _cpp: UnionCPPGreenTensor
    _cpp_type: ClassVar[UnionTypeCPPGreenTensor]

    def __init__(self) -> None:
        self._cpp = self._cpp_type()

    @overload
    def create_entries(
        self, kappa1: int, kappa2: int, tensor: "NDArray", tensor_unit: Optional[str] = None
    ) -> None: ...

    @overload
    def create_entries(
        self,
        kappa1: int,
        kappa2: int,
        tensor: Sequence["NDArray"],
        tensor_unit: Optional[str],
        omegas: Sequence[float],
        omega_unit: Optional[str] = None,
    ) -> None: ...

    @overload
    def create_entries(
        self,
        kappa1: int,
        kappa2: int,
        tensor: Sequence["NDArray"],
        *,
        omegas: Sequence[float],
        omega_unit: Optional[str] = None,
    ) -> None: ...

    def create_entries(
        self,
        kappa1: int,
        kappa2: int,
        tensor: Union["NDArray", Sequence["NDArray"]],
        tensor_unit: Optional[str] = None,
        omegas: Optional[Sequence[float]] = None,
        omega_unit: Optional[str] = None,
    ) -> None:
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
            omega_unit: The unit of the frequencies.
                Default None, which means that the frequencies must be given as pint object.


        """
        dimension: list[Dimension] = self.get_unit_dimension(kappa1, kappa2)
        if omegas is None:
            tensor_au = QuantityArray.convert_user_to_au(np.array(tensor), tensor_unit, dimension)
            if tensor_au.shape != (3**kappa1, 3**kappa2) or tensor_au.ndim != 2:
                raise ValueError("The tensor must be a 2D array of shape (3**kappa1, 3**kappa2).")
            self._cpp.create_entries(kappa1, kappa2, tensor_au)
            return

        tensors_au = [QuantityArray.convert_user_to_au(np.array(t), tensor_unit, dimension) for t in tensor]
        if not all(t.ndim == 2 for t in tensors_au):
            raise ValueError("The tensor must be a list of 2D arrays.")
        if not all(t.shape == (3**kappa1, 3**kappa2) for t in tensors_au):
            raise ValueError("The tensors must be of shape (3**kappa1, 3**kappa2).")

        omegas_au = QuantityArray.convert_user_to_au(np.array(omegas), omega_unit, "energy")
        self._cpp.create_entries(kappa1, kappa2, tensors_au, omegas_au.tolist())

    @overload
    def get(
        self,
        kappa1: int,
        kappa2: int,
        omega: Optional[float] = None,
        omega_unit: Optional[str] = None,
        unit: None = None,
    ) -> "PintArray": ...

    @overload
    def get(
        self, kappa1: int, kappa2: int, omega: Optional[float] = None, omega_unit: Optional[str] = None, *, unit: str
    ) -> "NDArray": ...

    def get(
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
        entries_cpp = self._cpp.get_entries(kappa1, kappa2)
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

        return QuantityArray.convert_au_to_user(tensor_au, self.get_unit_dimension(kappa1, kappa2), unit)

    def get_unit_dimension(self, kappa1: int, kappa2: int) -> list["Dimension"]:
        return ["green_tensor_00", *["inverse_distance" for _ in range(kappa1 + kappa2 + 1)]]  # type: ignore [list-item]


class GreenTensorReal(GreenTensor):
    _cpp: _backend.GreenTensorReal
    _cpp_type = _backend.GreenTensorReal


class GreenTensorComplex(GreenTensor):
    _cpp: _backend.GreenTensorComplex
    _cpp_type = _backend.GreenTensorComplex
