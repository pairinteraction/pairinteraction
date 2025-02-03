from collections.abc import Iterable
from typing import TYPE_CHECKING, Any, ClassVar, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction.backend import _backend
from pairinteraction.backend._wrapped.basis.Basis import BasisBase
from pairinteraction.backend._wrapped.get_functions import get_cpp_operator_type, get_cpp_parity
from pairinteraction.backend._wrapped.ket.KetPair import (
    KetPair,
    KetPairComplexDouble,
    KetPairComplexFloat,
    KetPairDouble,
    KetPairFloat,
)
from pairinteraction.backend._wrapped.OperatorType import OperatorType
from pairinteraction.backend._wrapped.Parity import Parity
from pairinteraction.units import QuantityAbstract, QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity
    from scipy.sparse import csr_matrix

    from pairinteraction.backend._wrapped.basis.BasisAtom import BasisAtom
    from pairinteraction.backend._wrapped.ket.KetAtom import KetAtom
    from pairinteraction.backend._wrapped.system.SystemAtom import SystemAtom
    from pairinteraction.units import Array

    KetAtomOrBasisAtom = TypeVar("KetAtomOrBasisAtom", KetAtom, BasisAtom, covariant=True)

KetPairType = TypeVar("KetPairType", bound=KetPair)
KetPairLike = Union[KetPair, Iterable["KetAtom"]]
BasisPairLike = Union["BasisPair", Iterable["BasisAtom"]]
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
        systems: Iterable["SystemAtom"],
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
    def get_amplitudes(self, ket_or_basis: KetPairLike) -> "np.ndarray[Any,Any]": ...

    @overload
    def get_amplitudes(self, ket_or_basis: BasisPairLike) -> "csr_matrix": ...

    def get_amplitudes(self, ket_or_basis: Union[KetPairLike, BasisPairLike]):
        ket_or_basis_cpp: list[Any]
        if not isinstance(ket_or_basis, Iterable):
            ket_or_basis_cpp = [ket_or_basis._cpp]
        else:
            ket_or_basis_cpp = [obj._cpp for obj in ket_or_basis]
        return self._cpp.get_amplitudes(*ket_or_basis_cpp)

    @overload
    def get_overlaps(self, ket_or_basis: KetPairLike) -> "np.ndarray[Any,Any]": ...

    @overload
    def get_overlaps(self, ket_or_basis: BasisPairLike) -> "csr_matrix": ...

    def get_overlaps(self, ket_or_basis: Union[KetPairLike, BasisPairLike]):
        ket_or_basis_cpp: list[Any]
        if not isinstance(ket_or_basis, Iterable):
            ket_or_basis_cpp = [ket_or_basis._cpp]
        else:
            ket_or_basis_cpp = [obj._cpp for obj in ket_or_basis]
        return self._cpp.get_overlaps(*ket_or_basis_cpp)

    @overload
    def get_matrix_elements(
        self, ket_or_basis: KetPairLike, operators: tuple[OperatorType, OperatorType], qs: tuple[int, int]
    ) -> "PlainQuantity[Array]": ...

    @overload
    def get_matrix_elements(
        self, ket_or_basis: KetPairLike, operators: tuple[OperatorType, OperatorType], qs: tuple[int, int], *, unit: str
    ) -> "Array": ...

    @overload
    def get_matrix_elements(
        self, ket_or_basis: BasisPairLike, operators: tuple[OperatorType, OperatorType], qs: tuple[int, int]
    ) -> "PlainQuantity[csr_matrix]": ...

    @overload
    def get_matrix_elements(
        self,
        ket_or_basis: BasisPairLike,
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        *,
        unit: str,
    ) -> "csr_matrix": ...

    def get_matrix_elements(
        self,
        ket_or_basis: Union[KetPairLike, BasisPairLike],
        operators: tuple[OperatorType, OperatorType],
        qs: tuple[int, int],
        *,
        unit: str = "pint",
    ):
        operators_cpp = [get_cpp_operator_type(operator) for operator in operators]
        ket_or_basis_cpp: list[Any]
        if not isinstance(ket_or_basis, Iterable):
            ket_or_basis_cpp = [ket_or_basis._cpp]
        else:
            ket_or_basis_cpp = [obj._cpp for obj in ket_or_basis]
        matrix_elements_au = self._cpp.get_matrix_elements(*ket_or_basis_cpp, *operators_cpp, *qs)
        matrix_elements: QuantityAbstract
        if isinstance(matrix_elements_au, np.ndarray):
            matrix_elements = QuantityArray.from_base(matrix_elements_au, operators)
        else:  # csr_matrix
            matrix_elements = QuantitySparse.from_base(matrix_elements_au, operators)
        return matrix_elements.to_unit(unit)


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
