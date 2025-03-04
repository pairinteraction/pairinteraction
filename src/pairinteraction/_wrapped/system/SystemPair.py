from collections.abc import Collection
from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction._wrapped.basis.BasisPair import BasisPairComplex, BasisPairReal
from pairinteraction._wrapped.system.System import SystemBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pint.facets.plain import PlainQuantity
    from typing_extensions import Self

BasisType = TypeVar("BasisType", "BasisPairReal", "BasisPairComplex")
UnionCPPSystemPair = Union[_backend.SystemPairReal, _backend.SystemPairComplex]
UnionTypeCPPSystemPair = Union[type[_backend.SystemPairReal], type[_backend.SystemPairComplex]]


class SystemPairBase(SystemBase[BasisType]):
    _cpp: UnionCPPSystemPair
    _cpp_type: ClassVar[UnionTypeCPPSystemPair]

    def __init__(self, basis: BasisType) -> None:
        """Create a system object for a pair of atoms.

        Use the given BasisPair object to create the system object.
        You can set the distance (vector) between the atoms afterwards via the corresponding methods.

        Examples:
            >>> import pairinteraction.real as pi
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
            SystemPairReal(BasisPairReal object with 140 states and 140 kets, is_diagonal=False)
            >>> pair_system = pair_system.diagonalize()
            >>> eigenvalues = pair_system.get_eigenvalues(unit="GHz")
            >>> print(f"{eigenvalues[0] - pair_energy:.5f}")
            -2.19974

        Args:
            basis: The :class:`pairinteraction.real.BasisPair` object that describes the basis of the system.

        """
        super().__init__(basis)
        self._distance_vector_au = [0, 0, np.inf]

    def set_order(self: "Self", order: int) -> "Self":
        self._cpp.set_order(order)
        return self

    def set_distance(
        self: "Self", distance: Union[float, "PlainQuantity[float]"], unit: Optional[str] = None
    ) -> "Self":
        return self.set_distance_vector([0, 0, distance], unit)

    def set_distance_vector(
        self: "Self",
        distance: Union["PlainQuantity[NDArray[Any]]", Collection[Union[float, "PlainQuantity[float]"]]],
        unit: Optional[str] = None,
    ) -> "Self":
        distance_au = [QuantityScalar.from_pint_or_unit(v, unit, "DISTANCE").to_base_unit() for v in distance]
        self._cpp.set_distance_vector(distance_au)
        self._distance_vector_au = distance_au
        return self

    @overload
    def get_distance_vector(self) -> list["PlainQuantity[float]"]: ...

    @overload
    def get_distance_vector(self, unit: str) -> list[float]: ...

    def get_distance_vector(self, unit: Optional[str] = None):  # type: ignore
        distance_vector = [QuantityScalar.from_base_unit(d, "DISTANCE") for d in self._distance_vector_au]
        return [d.to_pint_or_unit(unit) for d in distance_vector]

    @overload
    def get_distance(self) -> "PlainQuantity[float]": ...

    @overload
    def get_distance(self, unit: str) -> float: ...

    def get_distance(self, unit: Optional[str] = None):  # type: ignore
        distance = np.linalg.norm(self._distance_vector_au)
        return QuantityScalar.from_base_unit(float(distance), "DISTANCE").to_pint_or_unit(unit)


class SystemPairReal(SystemPairBase[BasisPairReal]):
    _cpp: _backend.SystemPairReal  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairReal
    _TypeBasis = BasisPairReal


class SystemPairComplex(SystemPairBase[BasisPairComplex]):
    _cpp: _backend.SystemPairComplex  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_type = _backend.SystemPairComplex
    _TypeBasis = BasisPairComplex


SystemPair = SystemPairBase[Any]
