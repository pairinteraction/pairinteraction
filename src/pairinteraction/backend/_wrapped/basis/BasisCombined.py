from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union, overload

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.Basis import BasisBase
from pairinteraction.backend._wrapped.ket.KetCombined import (
    KetCombinedBase,
    KetCombinedComplexDouble,
    KetCombinedComplexFloat,
    KetCombinedDouble,
    KetCombinedFloat,
)
from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtomBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    import scipy.sparse
    from numpy.typing import ArrayLike
    from pint.facets.plain import PlainQuantity

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtomBase
    from pairinteraction.backend._wrapped.ket.KetAtom import KetAtomBase

    KetAtomOrBasisAtom_t = TypeVar("KetAtomOrBasisAtom_t", "KetAtomBase", "BasisAtomBase[Any]", covariant=True)

Ket_t = TypeVar("Ket_t", bound=KetCombinedBase)
UnionCPPBasisCombined = Union[
    _backend.BasisCombinedFloat,
    _backend.BasisCombinedComplexFloat,
    _backend.BasisCombinedDouble,
    _backend.BasisCombinedComplexDouble,
]
UnionTypeCPPBasisCombinedCreator = Union[
    type[_backend.BasisCombinedCreatorFloat],
    type[_backend.BasisCombinedCreatorComplexFloat],
    type[_backend.BasisCombinedCreatorDouble],
    type[_backend.BasisCombinedCreatorComplexDouble],
]


class BasisCombinedBase(BasisBase[Ket_t]):
    _cpp: UnionCPPBasisCombined
    _cpp_creator: ClassVar[UnionTypeCPPBasisCombinedCreator]

    def __init__(
        self,
        systems: Sequence[SystemAtomBase[Any]],
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
        energy_unit: str = "pint",
    ) -> None:
        creator = self._cpp_creator()
        for system in systems:
            creator.add(system._cpp)  # type: ignore [reportPrivateUsage, arg-type]
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if energy is not None:
            min_energy_au = QuantityScalar(energy[0], energy_unit).to_base("energy")
            max_energy_au = QuantityScalar(energy[1], energy_unit).to_base("energy")
            creator.restrict_energy(min_energy_au, max_energy_au)
        self._cpp = creator.create()

    @overload
    def get_overlaps_from_product(
        self, ket_or_basis_1: "KetAtomBase", ket_or_basis_2: "KetAtomBase"
    ) -> "ArrayLike": ...

    @overload
    def get_overlaps_from_product(
        self, ket_or_basis_1: "BasisAtomBase[Any]", ket_or_basis_2: "BasisAtomBase[Any]"
    ) -> "scipy.sparse.csr_matrix": ...

    def get_overlaps_from_product(self, ket_or_basis_1: "KetAtomOrBasisAtom_t", ket_or_basis_2: "KetAtomOrBasisAtom_t"):  # type: ignore [reportUnknownParameterType]
        overlaps = self._cpp.get_overlaps(ket_or_basis_1._cpp, ket_or_basis_2._cpp)  # type: ignore
        return overlaps  # type: ignore [reportUnknownVariableType]


class BasisCombinedFloat(BasisCombinedBase[KetCombinedFloat]):
    _cpp: _backend.BasisCombinedFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisCombinedCreatorFloat
    _TypeKet = KetCombinedFloat


class BasisCombinedComplexFloat(BasisCombinedBase[KetCombinedComplexFloat]):
    _cpp: _backend.BasisCombinedComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisCombinedCreatorComplexFloat
    _TypeKet = KetCombinedComplexFloat


class BasisCombinedDouble(BasisCombinedBase[KetCombinedDouble]):
    _cpp: _backend.BasisCombinedDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisCombinedCreatorDouble
    _TypeKet = KetCombinedDouble


class BasisCombinedComplexDouble(BasisCombinedBase[KetCombinedComplexDouble]):
    _cpp: _backend.BasisCombinedComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisCombinedCreatorComplexDouble
    _TypeKet = KetCombinedComplexDouble
