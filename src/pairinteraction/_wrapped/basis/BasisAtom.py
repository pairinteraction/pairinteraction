from typing import TYPE_CHECKING, Any, ClassVar, Optional, Union, overload

import numpy as np

from pairinteraction import _backend
from pairinteraction._wrapped.basis.Basis import BasisBase
from pairinteraction._wrapped.cpp_types import OperatorType, Parity, get_cpp_operator_type, get_cpp_parity
from pairinteraction._wrapped.database.Database import Database
from pairinteraction._wrapped.ket.KetAtom import KetAtom
from pairinteraction.units import QuantityAbstract, QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction.units import Array

UnionCPPBasisAtom = Union[_backend.BasisAtomReal, _backend.BasisAtomComplex]
UnionTypeCPPBasisAtomCreator = Union[type[_backend.BasisAtomCreatorReal], type[_backend.BasisAtomCreatorComplex]]


class BasisAtomBase(BasisBase[KetAtom]):
    _cpp: UnionCPPBasisAtom
    _cpp_creator: ClassVar[UnionTypeCPPBasisAtomCreator]

    def __init__(  # noqa: C901
        self,
        species: str,
        n: Optional[tuple[int, int]] = None,
        nu: Optional[tuple[float, float]] = None,
        l: Optional[tuple[float, float]] = None,
        s: Optional[tuple[float, float]] = None,
        j: Optional[tuple[float, float]] = None,
        f: Optional[tuple[float, float]] = None,
        m: Optional[tuple[float, float]] = None,
        energy: Union[tuple[float, float], tuple["PlainQuantity[float]", "PlainQuantity[float]"], None] = None,
        energy_unit: Optional[str] = None,
        parity: Optional[Parity] = None,
        database: Optional[Database] = None,
        additional_kets: Optional[list[KetAtom]] = None,
    ) -> None:
        """Create a basis for a single atom.

        Add all KetAtom objects that match the given quantum numbers to the basis.
        The initial coefficients matrix is a unit matrix, i.e. the first basis state is the first ket, etc.
        The BasisAtom coefficients matrix will always be square,
        i.e. the number of kets is equal to the number of states.

        Examples:
            >>> import pairinteraction.real as pi
            >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
            >>> energy_min, energy_max = ket.get_energy(unit="GHz") - 100, ket.get_energy(unit="GHz") + 100
            >>> basis = pi.BasisAtom("Rb", n=(57, 63), l=(0, 3), energy=(energy_min, energy_max), energy_unit="GHz")
            >>> print(basis)
            BasisAtomReal object with 140 states and 140 kets

        Args:
            species: The species of the atom.
            n: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            nu: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            l: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            s: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            j: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            f: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            m: tuple of (min, max) values for this quantum number. Default None, i.e. add all available states.
            energy: tuple of (min, max) value for the energy. Default None, i.e. add all available states.
            energy_unit: In which unit the energy values are given, e.g. "GHz".
                Default None, i.e. energy is provided as pint object.
            parity: The parity of the states to consider. Default None, i.e. add all available states.
            database: Which database to use. Default None, i.e. use the global database instance.
            additional_kets: List of additional kets to add to the basis. Default None.

        """
        creator = self._cpp_creator()
        creator.set_species(species)
        if f is not None:
            creator.restrict_quantum_number_f(*f)
        if m is not None:
            creator.restrict_quantum_number_m(*m)
        if parity is not None:
            creator.restrict_parity(get_cpp_parity(parity))
        if n is not None:
            creator.restrict_quantum_number_n(*n)
        if nu is not None:
            creator.restrict_quantum_number_nu(*nu)
        if l is not None:
            creator.restrict_quantum_number_l(*l)
        if s is not None:
            creator.restrict_quantum_number_s(*s)
        if j is not None:
            creator.restrict_quantum_number_j(*j)
        if energy is not None:
            min_energy_au = QuantityScalar.from_pint_or_unit(energy[0], energy_unit, "ENERGY").to_base_unit()
            max_energy_au = QuantityScalar.from_pint_or_unit(energy[1], energy_unit, "ENERGY").to_base_unit()
            creator.restrict_energy(min_energy_au, max_energy_au)
        if database is None:
            if Database.get_global_database() is None:
                Database.initialize_global_database()
            database = Database.get_global_database()
        if additional_kets is not None:
            for ket in additional_kets:
                creator.append_ket(ket._cpp)  # type: ignore [reportPrivateUsage]
        self._cpp = creator.create(database._cpp)  # type: ignore [reportPrivateUsage]

    @overload
    def get_amplitudes(self, ket_or_basis: KetAtom) -> "np.ndarray[Any,Any]": ...

    @overload
    def get_amplitudes(self, ket_or_basis: "Self") -> "csr_matrix": ...

    def get_amplitudes(self, ket_or_basis: Union[KetAtom, "Self"]):
        return self._cpp.get_amplitudes(ket_or_basis._cpp)

    @overload
    def get_overlaps(self, ket_or_basis: KetAtom) -> "np.ndarray[Any,Any]": ...

    @overload
    def get_overlaps(self, ket_or_basis: "Self") -> "csr_matrix": ...

    def get_overlaps(self, ket_or_basis: Union[KetAtom, "Self"]):
        return self._cpp.get_overlaps(ket_or_basis._cpp)

    @overload
    def get_matrix_elements(self, ket_or_basis: KetAtom, operator: OperatorType, q: int) -> "PlainQuantity[Array]": ...

    @overload
    def get_matrix_elements(self, ket_or_basis: KetAtom, operator: OperatorType, q: int, unit: str) -> "Array": ...

    @overload
    def get_matrix_elements(
        self, ket_or_basis: "Self", operator: OperatorType, q: int
    ) -> "PlainQuantity[csr_matrix]": ...

    @overload
    def get_matrix_elements(self, ket_or_basis: "Self", operator: OperatorType, q: int, unit: str) -> "csr_matrix": ...

    def get_matrix_elements(
        self, ket_or_basis: Union[KetAtom, "Self"], operator: OperatorType, q: int, unit: Optional[str] = None
    ):
        matrix_elements_au = self._cpp.get_matrix_elements(ket_or_basis._cpp, get_cpp_operator_type(operator), q)
        matrix_elements: QuantityAbstract
        if isinstance(matrix_elements_au, np.ndarray):
            matrix_elements = QuantityArray.from_base_unit(matrix_elements_au, operator)
        else:  # csr_matrix
            matrix_elements = QuantitySparse.from_base_unit(matrix_elements_au, operator)
        return matrix_elements.to_pint_or_unit(unit)


class BasisAtomReal(BasisAtomBase):
    _cpp: _backend.BasisAtomReal  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorReal
    _TypeKet = KetAtom


class BasisAtomComplex(BasisAtomBase):
    _cpp: _backend.BasisAtomComplex  # type: ignore [reportIncompatibleVariableOverride]
    _cpp_creator = _backend.BasisAtomCreatorComplex
    _TypeKet = KetAtom


BasisAtom = BasisAtomBase
