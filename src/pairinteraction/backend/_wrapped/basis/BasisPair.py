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

KetPairType = TypeVar("KetPairType", bound=KetPair)
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


class BasisPairBase(BasisBase[KetPairType]):
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
        """Create a basis for a pair of atoms.

        Add all product states of the eigenstates of two given SystemAtom objects to the basis,
        which pair energy is within the given energy range.
        You can also specify which total magnetic quantum number m the pair should have (if it is conserved)
        and the product of the parities of the two atoms.
        Due to the possible restrictions of the basis states, the BasisPair coefficients matrix will in general
        not be square but (n x d),
        where n is the number of all involved kets (typically basis1.number_of_kets * basis2.number_of_kets)
        and d is the number of basis states (after applying the restrictions).

        Examples:
            >>> import pairinteraction.backend.double as pi
            >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
            >>> system = pi.SystemAtom(basis).set_electric_field([0.1, 0, 0.1], unit="V/cm").diagonalize()
            >>> pair_energy = 2 * pi.calculate_energy(ket, system, unit="GHz")
            >>> pair_basis = pi.BasisPair(
            ...     [system, system],
            ...     energy=(pair_energy - 3, pair_energy + 3),
            ...     energy_unit="GHz",
            ... )
            >>> print(pair_basis)
            BasisPairDouble object with 140 states and 140 kets

        Args:
            systems: tuple of two SystemAtom objects, which define the two atoms, from which the BasisPair is build.
                Both systems have to be diagonalized before creating the BasisPair.
            m: tuple of (min, max) values for the total magnetic quantum number m of the pair state.
                Default None, i.e. no restriction.
            product_of_parities: The product parity of the states to consider.
                Default None, i.e. add all available states.
            energy: tuple of (min, max) value for the pair energy. Default None, i.e. add all available states.
            energy_unit: In which unit the energy values are given, e.g. "GHz".
                Default "pint", i.e. energy is provided as pint object.

        """
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
