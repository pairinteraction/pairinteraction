import logging
from collections.abc import Collection
from typing import TYPE_CHECKING, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction._wrapped.basis.basis_atom import BasisAtomComplex, BasisAtomReal
from pairinteraction._wrapped.system.system import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.basis_atom import BasisAtom
    from pairinteraction._wrapped.ket.ket_atom import KetAtom
    from pairinteraction.units import PintArray, PintFloat

BasisType = TypeVar("BasisType", bound="BasisAtom", covariant=True)
UnionCPPSystemAtom = Union[_backend.SystemAtomReal, _backend.SystemAtomComplex]
UnionTypeCPPSystemAtom = Union[type[_backend.SystemAtomReal], type[_backend.SystemAtomComplex]]

logger = logging.getLogger(__name__)


class SystemAtom(SystemBase[BasisType]):
    """System of a single atom.

    Use the given BasisAtom object to create the system object.
    You can set the electric and magnetic fields and enable diamagnetism afterwards via the corresponding methods.

    Examples:
        >>> import pairinteraction.real as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> system = pi.SystemAtom(basis)
        >>> system = system.set_magnetic_field([0, 0, 1], unit="gauss")
        >>> system = system.set_electric_field([0.1, 0, 0.1], unit="V/cm")
        >>> system = system.enable_diamagnetism(True)
        >>> print(system)
        SystemAtom(BasisAtom(|Rb:58,S_1/2,-1/2⟩ ... |Rb:63,F_5/2,5/2⟩), is_diagonal=False)
        >>> system = system.diagonalize()
        >>> eigenenergies = system.get_eigenenergies(unit="GHz")
        >>> print(f"{eigenenergies[0] - ket.get_energy(unit='GHz'):.5f}")
        -75.51823

    """

    _cpp: UnionCPPSystemAtom
    _cpp_type: ClassVar[UnionTypeCPPSystemAtom]

    def __init__(self, basis: BasisType) -> None:
        """Create a system object for a single atom.

        Args:
            basis: The :class:`pairinteraction.real.BasisAtom` object that describes the basis of the system.

        """
        super().__init__(basis)

    def set_electric_field(
        self: "Self",
        electric_field: Union["PintArray", Collection[Union[float, "PintFloat"]]],
        unit: Optional[str] = None,
    ) -> "Self":
        electric_field_au = [
            QuantityScalar.from_pint_or_unit(v, unit, "electric_field").to_base_unit() for v in electric_field
        ]
        self._cpp.set_electric_field(electric_field_au)
        return self

    def set_magnetic_field(
        self: "Self",
        magnetic_field: Union["PintArray", Collection[Union[float, "PintFloat"]]],
        unit: Optional[str] = None,
    ) -> "Self":
        magnetic_field_au = [
            QuantityScalar.from_pint_or_unit(v, unit, "magnetic_field").to_base_unit() for v in magnetic_field
        ]
        self._cpp.set_magnetic_field(magnetic_field_au)
        return self

    def enable_diamagnetism(self: "Self", enable: bool = True) -> "Self":
        self._cpp.enable_diamagnetism(enable)
        return self

    @overload
    def get_corresponding_energy(self: "Self", ket: "KetAtom", unit: None = None) -> "PintFloat": ...

    @overload
    def get_corresponding_energy(self: "Self", ket: "KetAtom", unit: str) -> float: ...

    def get_corresponding_energy(self: "Self", ket: "KetAtom", unit: Optional[str] = None) -> Union[float, "PintFloat"]:
        overlaps = self.get_eigenbasis().get_overlaps(ket)
        idx = np.argmax(overlaps)
        if overlaps[idx] <= 0.5:
            logger.warning(
                "The provided ket states does not correspond to an eigenstate of the system in a unique way."
            )
        return self.get_eigenenergies(unit=unit)[idx]  # type: ignore [index,no-any-return] # PintArray does not know it can be indexed


class SystemAtomReal(SystemAtom[BasisAtomReal]):
    _cpp: _backend.SystemAtomReal
    _cpp_type = _backend.SystemAtomReal
    _TypeBasis = BasisAtomReal


class SystemAtomComplex(SystemAtom[BasisAtomComplex]):
    _cpp: _backend.SystemAtomComplex
    _cpp_type = _backend.SystemAtomComplex
    _TypeBasis = BasisAtomComplex
