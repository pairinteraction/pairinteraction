# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Sequence
from typing import TYPE_CHECKING, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction import _backend

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.units import NDArray, PintFloat

    Quantity = TypeVar("Quantity", float, "PintFloat")

UnionCPPGreenTensor = Union[_backend.GreenTensorReal, _backend.GreenTensorComplex]
UnionTypeCPPGreenTensor = Union[type[_backend.GreenTensorReal], type[_backend.GreenTensorComplex]]

UnionCPPConstantEntry = Union[_backend.ConstantEntryReal, _backend.ConstantEntryComplex]
UnionCPPOmegaDependentEntry = Union[_backend.OmegaDependentEntryReal, _backend.OmegaDependentEntryComplex]


class GreenTensor:
    _cpp: UnionCPPGreenTensor
    _cpp_type: ClassVar[UnionTypeCPPGreenTensor]

    def __init__(self) -> None:
        self._cpp = self._cpp_type()

    @overload
    def set_entries(self, kappa1: int, kappa2: int, tensor: "NDArray") -> None: ...

    @overload
    def set_entries(self, kappa1: int, kappa2: int, tensor: Sequence["NDArray"], omegas: Sequence[float]) -> None: ...

    def set_entries(
        self,
        kappa1: int,
        kappa2: int,
        tensor: Union["NDArray", Sequence["NDArray"]],
        omegas: Optional[Sequence[float]] = None,
    ) -> None:
        """Set the entries of the Green tensor.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.
            tensor: The tensor to set.
                Either a single tensor to set a constant green tensor
                or a list of tensors to set a frequency-dependent green tensor.
            omegas: Only needed if a list of tensors is given.
                The frequencies of the tensors.

        """
        if omegas is None:
            tensor = np.array(tensor)
            if tensor.ndim != 2:
                raise ValueError("The tensor must be a 2D array.")
            self._cpp.set_entries(kappa1, kappa2, tensor)
        else:
            tensor = [np.array(t) for t in tensor]
            if any(t.ndim != 2 for t in tensor):
                raise ValueError("The tensor must be a list of 2D arrays.")
            self._cpp.set_entries(kappa1, kappa2, tensor, omegas)

    def get_entries(self, kappa1: int, kappa2: int) -> list[Union["ConstantEntry", "OmegaDependentEntry"]]:
        """Get the entries of the Green tensor for the given indices kappa1 and kappa2.

        Args:
            kappa1: The rank of the first multipole operator.
            kappa2: The rank of the second multipole operator.

        Returns:
            A list of entries of the Green tensor.

        """
        entries_cpp = self._cpp.get_entries(kappa1, kappa2)
        entries: list[Union[ConstantEntry, OmegaDependentEntry]] = []
        for entry_cpp in entries_cpp:
            if isinstance(entry_cpp, _backend.ConstantEntryReal):
                entries.append(ConstantEntryReal._from_cpp_object(entry_cpp))
            elif isinstance(entry_cpp, _backend.ConstantEntryComplex):
                entries.append(ConstantEntryComplex._from_cpp_object(entry_cpp))
            elif isinstance(entry_cpp, _backend.OmegaDependentEntryReal):
                entries.append(OmegaDependentEntryReal._from_cpp_object(entry_cpp))
            elif isinstance(entry_cpp, _backend.OmegaDependentEntryComplex):
                entries.append(OmegaDependentEntryComplex._from_cpp_object(entry_cpp))
            else:
                raise TypeError(f"Unknown entry type: {type(entry_cpp)}")
        return entries


class GreenTensorReal(GreenTensor):
    _cpp: _backend.GreenTensorReal
    _cpp_type = _backend.GreenTensorReal


class GreenTensorComplex(GreenTensor):
    _cpp: _backend.GreenTensorComplex
    _cpp_type = _backend.GreenTensorComplex


class ConstantEntry:
    _cpp: UnionCPPConstantEntry

    def __init__(self) -> None:
        raise NotImplementedError("ConstantEntry objects cannot be created directly.")

    @classmethod
    def _from_cpp_object(cls: "type[Self]", cpp_obj: UnionCPPConstantEntry) -> "Self":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def row(self) -> int:
        """Return the row index of the entry."""
        return self._cpp.row()

    def col(self) -> int:
        """Return the column index of the entry."""
        return self._cpp.col()

    def val(self) -> complex:
        """Return the constant value of the entry."""
        return self._cpp.val()


class ConstantEntryReal(ConstantEntry):
    _cpp: _backend.ConstantEntryReal


class ConstantEntryComplex(ConstantEntry):
    _cpp: _backend.ConstantEntryComplex


class OmegaDependentEntry:
    _cpp: UnionCPPOmegaDependentEntry

    def __init__(self) -> None:
        raise NotImplementedError("ConstantEntry objects cannot be created directly.")

    @classmethod
    def _from_cpp_object(cls: "type[Self]", cpp_obj: UnionCPPOmegaDependentEntry) -> "Self":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        return obj

    def row(self) -> int:
        """Return the row index of the entry."""
        return self._cpp.row()

    def col(self) -> int:
        """Return the column index of the entry."""
        return self._cpp.col()

    def val(self, omega: float) -> Union[float, complex]:
        """Return the omega dependent value of the entry.

        Args:
            omega: The frequency at which to evaluate the green tensor entry.

        """
        return self._cpp.val(omega)


class OmegaDependentEntryReal(OmegaDependentEntry):
    _cpp: _backend.OmegaDependentEntryReal


class OmegaDependentEntryComplex(OmegaDependentEntry):
    _cpp: _backend.OmegaDependentEntryComplex
