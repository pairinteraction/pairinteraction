import logging
from collections.abc import Collection
from typing import TYPE_CHECKING, Any, ClassVar, TypeVar, Union, overload

import numpy as np

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.BasisAtom import (
    BasisAtomComplexDouble,
    BasisAtomComplexFloat,
    BasisAtomDouble,
    BasisAtomFloat,
)
from pairinteraction.backend._wrapped.system.System import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity
    from typing_extensions import Self

    from pairinteraction.backend._wrapped.ket.KetAtom import KetAtom
    from pairinteraction.units import Array

BasisType = TypeVar("BasisType", "BasisAtomFloat", "BasisAtomComplexFloat", "BasisAtomDouble", "BasisAtomComplexDouble")
UnionCPPSystemAtom = Union[
    _backend.SystemAtomFloat,
    _backend.SystemAtomComplexFloat,
    _backend.SystemAtomDouble,
    _backend.SystemAtomComplexDouble,
]
UnionTypeCPPSystemAtom = Union[
    type[_backend.SystemAtomFloat],
    type[_backend.SystemAtomComplexFloat],
    type[_backend.SystemAtomDouble],
    type[_backend.SystemAtomComplexDouble],
]

logger = logging.getLogger(__name__)


class SystemAtomBase(SystemBase[BasisType]):
    _cpp: UnionCPPSystemAtom
    _cpp_type: ClassVar[UnionTypeCPPSystemAtom]

    def __init__(self, basis: BasisType) -> None:
        """Create a system object for a single atom.

        Use the given BasisAtom object to create the system object.
        You can set the electric and magnetic fields and enable diamagnetism afterwards via the corresponding methods.

        Examples:
            >>> import pairinteraction.backend.double as pi
            >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
            >>> system = pi.SystemAtom(basis)
            >>> system = system.set_magnetic_field([0, 0, 1], unit="gauss")
            >>> system = system.set_electric_field([0.1, 0, 0.1], unit="V/cm")
            >>> system = system.enable_diamagnetism(True)
            >>> print(system)
            SystemAtomDouble(BasisAtomDouble object with 192 states and 192 kets, is_diagonal=False)
            >>> system = system.diagonalize()
            >>> eigenvalues = system.get_eigenvalues(unit="GHz")
            >>> print(f"{eigenvalues[0] - ket.get_energy(unit='GHz'):.5f}")
            -75.51823

        Args:
            basis: The :class:`pairinteraction.backend.double.BasisAtom` object that describes the basis of the system.

        """
        super().__init__(basis)

    def set_electric_field(
        self: "Self",
        electric_field: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ) -> "Self":
        electric_field_au = [QuantityScalar(v, unit).to_base("ELECTRIC_FIELD") for v in electric_field]
        self._cpp.set_electric_field(electric_field_au)
        return self

    def set_magnetic_field(
        self: "Self",
        magnetic_field: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: str = "pint",
    ) -> "Self":
        magnetic_field_au = [QuantityScalar(v, unit).to_base("MAGNETIC_FIELD") for v in magnetic_field]
        self._cpp.set_magnetic_field(magnetic_field_au)
        return self

    def enable_diamagnetism(self: "Self", enable: bool = True) -> "Self":
        self._cpp.enable_diamagnetism(enable)
        return self

    @overload
    def get_corresponding_energy(self: "Self", ket: "KetAtom") -> "PlainQuantity[float]": ...

    @overload
    def get_corresponding_energy(self: "Self", ket: "KetAtom", *, unit: str) -> float: ...

    def get_corresponding_energy(self: "Self", ket: "KetAtom", *, unit: str = "pint"):
        overlaps = self.get_eigenbasis().get_overlaps(ket)
        idx = np.argmax(overlaps)
        if overlaps[idx] <= 0.5:
            logger.warning(
                "The provided ket states does not correspond to an eigenstate of the system in a unique way."
            )
        return self.get_eigenvalues(unit=unit)[idx]


class SystemAtomFloat(SystemAtomBase[BasisAtomFloat]):
    _cpp: _backend.SystemAtomFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomFloat
    _TypeBasis = BasisAtomFloat


class SystemAtomComplexFloat(SystemAtomBase[BasisAtomComplexFloat]):
    _cpp: _backend.SystemAtomComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomComplexFloat
    _TypeBasis = BasisAtomComplexFloat


class SystemAtomDouble(SystemAtomBase[BasisAtomDouble]):
    _cpp: _backend.SystemAtomDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomDouble
    _TypeBasis = BasisAtomDouble


class SystemAtomComplexDouble(SystemAtomBase[BasisAtomComplexDouble]):
    _cpp: _backend.SystemAtomComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemAtomComplexDouble
    _TypeBasis = BasisAtomComplexDouble


SystemAtom = SystemAtomBase[Any]
