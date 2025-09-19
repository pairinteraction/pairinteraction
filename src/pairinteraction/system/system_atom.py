# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction.basis import BasisAtom, BasisAtomReal
from pairinteraction.system.system_base import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.ket import KetAtom
    from pairinteraction.units import (
        ArrayLike,
        PintArray,  # noqa: F401  # needed for sphinx to recognize PintArrayLike
        PintArrayLike,
        PintFloat,
    )


logger = logging.getLogger(__name__)


class SystemAtom(SystemBase[BasisAtom]):
    """System of a single atom.

    Use the given BasisAtom object to create the system object.
    You can set the electric and magnetic fields and enable diamagnetism afterwards via the corresponding methods.

    Examples:
        >>> import pairinteraction as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> system = pi.SystemAtom(basis)
        >>> system.set_magnetic_field([0, 0, 1], unit="gauss").set_electric_field([0.1, 0, 0.1], unit="V/cm")
        SystemAtom(BasisAtom(n=(58, 63), l=(0, 3)), is_diagonal=False)
        >>> system.diagonalize()
        SystemAtom(BasisAtom(|Rb:58,S_1/2,-1/2⟩ ... |Rb:63,F_5/2,5/2⟩), is_diagonal=True)
        >>> eigenenergies = system.get_eigenenergies(unit="GHz")
        >>> print(f"{eigenenergies[0] - ket.get_energy(unit='GHz'):.5f}")
        -75.51823

    """

    _cpp: _backend.SystemAtomComplex
    _cpp_type = _backend.SystemAtomComplex
    _basis_class = BasisAtom

    def __init__(self, basis: BasisAtom) -> None:
        """Create a system object for a single atom.

        Args:
            basis: The :class:`pairinteraction.BasisAtom` object that describes the basis of the system.

        """
        self._cpp = self._cpp_type(basis._cpp)
        self._basis = basis

    def set_electric_field(
        self: Self,
        electric_field: PintArrayLike | ArrayLike,
        unit: str | None = None,
    ) -> Self:
        """Set the electric field for the system.

        Args:
            electric_field: The electric field to set for the system.
            unit: The unit of the electric field, e.g. "V/cm".
                Default None expects a `pint.Quantity`.

        """
        electric_field_au = [QuantityScalar.convert_user_to_au(v, unit, "electric_field") for v in electric_field]
        self._cpp.set_electric_field(electric_field_au)
        return self

    def set_magnetic_field(
        self: Self,
        magnetic_field: PintArrayLike | ArrayLike,
        unit: str | None = None,
    ) -> Self:
        """Set the magnetic field for the system.

        Args:
            magnetic_field: The magnetic field to set for the system.
            unit: The unit of the magnetic field, e.g. "gauss".
                Default None expects a `pint.Quantity`.

        """
        magnetic_field_au = [QuantityScalar.convert_user_to_au(v, unit, "magnetic_field") for v in magnetic_field]
        self._cpp.set_magnetic_field(magnetic_field_au)
        return self

    def set_diamagnetism_enabled(self: Self, enable: bool = True) -> Self:
        """Enable or disable diamagnetism for the system.

        Args:
            enable: Whether to enable or disable diamagnetism.

        """
        self._cpp.set_diamagnetism_enabled(enable)
        return self

    def set_distance_to_ion(
        self: Self, distance: float | PintFloat, angle_degree: float = 0, unit: str | None = None
    ) -> Self:
        distance_vector = [np.sin(np.deg2rad(angle_degree)) * distance, 0, np.cos(np.deg2rad(angle_degree)) * distance]
        return self.set_ion_distance_vector(distance_vector, unit)

    def set_ion_distance_vector(
        self: Self,
        distance: PintArrayLike | ArrayLike,
        unit: str | None = None,
    ) -> Self:
        distance_au = [QuantityScalar.convert_user_to_au(v, unit, "distance") for v in distance]
        self._cpp.set_ion_distance_vector(distance_au)
        return self

    def set_ion_charge(self: Self, charge: float | PintFloat, unit: str | None = None) -> Self:
        charge_au = QuantityScalar.convert_user_to_au(charge, unit, "charge")
        self._cpp.set_ion_charge(charge_au)
        return self

    def set_ion_interaction_order(self: Self, order: int) -> Self:
        self._cpp.set_ion_interaction_order(order)
        return self

    @overload
    def get_corresponding_energy(self: Self, ket: KetAtom, unit: None = None) -> PintFloat: ...

    @overload
    def get_corresponding_energy(self: Self, ket: KetAtom, unit: str) -> float: ...

    def get_corresponding_energy(self: Self, ket: KetAtom, unit: str | None = None) -> float | PintFloat:
        overlaps = self.get_eigenbasis().get_overlaps(ket)
        idx = np.argmax(overlaps)
        if overlaps[idx] <= 0.5:
            logger.warning(
                "The provided ket states does not correspond to an eigenstate of the system in a unique way."
            )
        return self.get_eigenenergies(unit=unit)[idx]  # type: ignore [index,no-any-return] # PintArray does not know it can be indexed


class SystemAtomReal(SystemAtom):
    _cpp: _backend.SystemAtomReal  # type: ignore [assignment]
    _cpp_type = _backend.SystemAtomReal  # type: ignore [assignment]
    _basis_class = BasisAtomReal
