from collections.abc import Collection
from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.BasisPair import (
    BasisPairComplexDouble,
    BasisPairComplexFloat,
    BasisPairDouble,
    BasisPairFloat,
)
from pairinteraction.backend._wrapped.system.System import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity
    from typing_extensions import Self

    from pairinteraction.units import Array

BasisType = TypeVar("BasisType", "BasisPairFloat", "BasisPairComplexFloat", "BasisPairDouble", "BasisPairComplexDouble")
UnionCPPSystemPair = Union[
    _backend.SystemPairFloat,
    _backend.SystemPairComplexFloat,
    _backend.SystemPairDouble,
    _backend.SystemPairComplexDouble,
]
UnionTypeCPPSystemPair = Union[
    type[_backend.SystemPairFloat],
    type[_backend.SystemPairComplexFloat],
    type[_backend.SystemPairDouble],
    type[_backend.SystemPairComplexDouble],
]


class SystemPairBase(SystemBase[BasisType]):
    _cpp: UnionCPPSystemPair
    _cpp_type: ClassVar[UnionTypeCPPSystemPair]

    def __init__(self, basis: BasisType) -> None:
        """Create a system object for a pair of atoms.

        Use the given BasisPair object to create the system object.
        You can set the distance (vector) between the atoms afterwards via the corresponding methods.

        Examples:
            >>> import pairinteraction.backend.double as pi
            >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
            >>> system = pi.SystemAtom(basis).set_electric_field([0.1, 0, 0.1], unit="V/cm").diagonalize()
            >>> pair_energy = 2 * system.get_corresponding_energy(ket, unit="GHz")
            >>> pair_basis = pi.BasisPair(
            ...     [system, system],
            ...     energy=(pair_energy - 3, pair_energy + 3),
            ...     energy_unit="GHz",
            ... )
            >>> pair_system = pi.SystemPair(pair_basis).set_distance(5, unit="micrometer").set_order(3)
            >>> print(pair_system)
            SystemPairDouble(BasisPairDouble object with 140 states and 140 kets, is_diagonal=False)
            >>> pair_system = pair_system.diagonalize()
            >>> eigenvalues = pair_system.get_eigenvalues(unit="GHz")
            >>> print(f"{eigenvalues[0] - pair_energy:.5f}")
            -2.19974

        Args:
            basis: The :class:`pairinteraction.backend.double.BasisPair` object that describes the basis of the system.

        """
        super().__init__(basis)

    def set_order(self: "Self", order: int) -> "Self":
        self._cpp.set_order(order)
        return self

    def set_distance(
        self: "Self", distance: Union[float, "PlainQuantity[float]"], unit: Optional[str] = None
    ) -> "Self":
        distance_au = QuantityScalar(distance, unit).to_base("DISTANCE")
        self._cpp.set_distance(distance_au)
        return self

    def set_distance_vector(
        self: "Self",
        distance: Union["PlainQuantity[Array]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: Optional[str] = None,
    ) -> "Self":
        distance_au = [QuantityScalar(v, unit).to_base("DISTANCE") for v in distance]
        self._cpp.set_distance_vector(distance_au)
        return self


class SystemPairFloat(SystemPairBase[BasisPairFloat]):
    _cpp: _backend.SystemPairFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairFloat
    _TypeBasis = BasisPairFloat


class SystemPairComplexFloat(SystemPairBase[BasisPairComplexFloat]):
    _cpp: _backend.SystemPairComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairComplexFloat
    _TypeBasis = BasisPairComplexFloat


class SystemPairDouble(SystemPairBase[BasisPairDouble]):
    _cpp: _backend.SystemPairDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairDouble
    _TypeBasis = BasisPairDouble


class SystemPairComplexDouble(SystemPairBase[BasisPairComplexDouble]):
    _cpp: _backend.SystemPairComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairComplexDouble
    _TypeBasis = BasisPairComplexDouble


SystemPair = SystemPairBase[Any]
