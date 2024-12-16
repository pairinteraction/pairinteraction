from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.Basis import BasisBase
from pairinteraction.backend._wrapped.ket.KetPair import (
    KetPair,
    KetPairComplexDouble,
    KetPairComplexFloat,
    KetPairDouble,
    KetPairFloat,
)
from pairinteraction.backend._wrapped.Parity import Parity, get_cpp_parity
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    import scipy.sparse
    from pint.facets.plain import PlainQuantity

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtom
    from pairinteraction.backend._wrapped.ket.KetAtom import KetAtom
    from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtom

    KetAtomOrBasisAtom = TypeVar("KetAtomOrBasisAtom", KetAtom, BasisAtom, covariant=True)

Ket_t = TypeVar("Ket_t", bound=KetPair)
UnionCPPBasisPair = Union[
    _backend.BasisPairFloat,
    _backend.BasisPairComplexFloat,
    _backend.BasisPairDouble,
    _backend.BasisPairComplexDouble,
]
UnionTypeCPPBasisPairCreator = Union[
    type[_backend.BasisPairCreatorFloat],
    type[_backend.BasisPairCreatorComplexFloat],
    type[_backend.BasisPairCreatorDouble],
    type[_backend.BasisPairCreatorComplexDouble],
]


class BasisPairBase(BasisBase[Ket_t]):
    _cpp: UnionCPPBasisPair
    _cpp_creator: ClassVar[UnionTypeCPPBasisPairCreator]

    def __init__(
        self,
        systems: Sequence["SystemAtom"],
        m: Optional[tuple[float, float]] = None,
        product_of_parities: Optional[Parity] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
        energy_unit: str = "pint",
    ) -> None:
        creator = self._cpp_creator()
        for system in systems:
            creator.add(system._cpp)  # type: ignore [reportPrivateUsage, arg-type]
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if product_of_parities is not None:
            creator.restrict_product_of_parities(get_cpp_parity(product_of_parities))
        if energy is not None:
            min_energy_au = QuantityScalar(energy[0], energy_unit).to_base("ENERGY")
            max_energy_au = QuantityScalar(energy[1], energy_unit).to_base("ENERGY")
            creator.restrict_energy(min_energy_au, max_energy_au)
        self._cpp = creator.create()

    @overload
    def get_amplitudes_with_product_state(
        self, ket_or_basis_1: "KetAtom", ket_or_basis_2: "KetAtom"
    ) -> "np.ndarray[Any, Any]": ...

    @overload
    def get_amplitudes_with_product_state(
        self, ket_or_basis_1: "BasisAtom", ket_or_basis_2: "BasisAtom"
    ) -> "scipy.sparse.csr_matrix": ...

    def get_amplitudes_with_product_state(
        self, ket_or_basis_1: "KetAtomOrBasisAtom", ket_or_basis_2: "KetAtomOrBasisAtom"
    ):  # type: ignore [reportUnknownParameterType]
        overlaps = self._cpp.get_amplitudes(ket_or_basis_1._cpp, ket_or_basis_2._cpp)  # type: ignore
        return overlaps  # type: ignore [reportUnknownVariableType]

    @overload
    def get_overlaps_with_product_state(
        self, ket_or_basis_1: "KetAtom", ket_or_basis_2: "KetAtom"
    ) -> "np.ndarray[Any, Any]": ...

    @overload
    def get_overlaps_with_product_state(
        self, ket_or_basis_1: "BasisAtom", ket_or_basis_2: "BasisAtom"
    ) -> "scipy.sparse.csr_matrix": ...

    def get_overlaps_with_product_state(
        self, ket_or_basis_1: "KetAtomOrBasisAtom", ket_or_basis_2: "KetAtomOrBasisAtom"
    ):  # type: ignore [reportUnknownParameterType]
        overlaps = self._cpp.get_overlaps(ket_or_basis_1._cpp, ket_or_basis_2._cpp)  # type: ignore
        return overlaps  # type: ignore [reportUnknownVariableType]


class BasisPairFloat(BasisPairBase[KetPairFloat]):
    _cpp: _backend.BasisPairFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisPairCreatorFloat
    _TypeKet = KetPairFloat


class BasisPairComplexFloat(BasisPairBase[KetPairComplexFloat]):
    _cpp: _backend.BasisPairComplexFloat  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisPairCreatorComplexFloat
    _TypeKet = KetPairComplexFloat


class BasisPairDouble(BasisPairBase[KetPairDouble]):
    _cpp: _backend.BasisPairDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisPairCreatorDouble
    _TypeKet = KetPairDouble


class BasisPairComplexDouble(BasisPairBase[KetPairComplexDouble]):
    _cpp: _backend.BasisPairComplexDouble  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisPairCreatorComplexDouble
    _TypeKet = KetPairComplexDouble


BasisPair = BasisPairBase[Any]
