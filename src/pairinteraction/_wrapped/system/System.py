from abc import ABC
from typing import TYPE_CHECKING, Any, ClassVar, Generic, Optional, TypeVar, Union, overload

from pairinteraction import _backend
from pairinteraction._wrapped.cpp_types import Diagonalizer, FloatType, get_cpp_diagonalizer
from pairinteraction.units import QuantityArray, QuantityScalar, QuantitySparse

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pint.facets.plain import PlainQuantity
    from scipy.sparse import csr_matrix
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.BasisAtom import BasisAtom
    from pairinteraction._wrapped.basis.BasisPair import BasisPair

    Quantity = TypeVar("Quantity", float, PlainQuantity[float])


BasisType = TypeVar("BasisType", bound=Union["BasisAtom", "BasisPair"])
UnionCPPSystem = Any
# UnionCPPSystem is supposed to be System(|System)(Atom|Pair)(Real|Complex)
UnionTypeCPPSystem = Any


class SystemBase(ABC, Generic[BasisType]):
    """Base class for all System objects.

    The system objects are meant to represent the full physical system, including the basis
    as well as possible fields and interactions which form the Hamiltonian.

    Typically, the system objects are simply constructed from a basis object.
    The fields and interactions can then be set afterwards.

    All system objects share a few common attributes and methods, that are defined in this base class, e.g.:
    - the basis of the system,
    - the matrix stored as scipy sparse matrix, which represents the Hamiltonian,
    - ...
    """

    _cpp: UnionCPPSystem
    _cpp_type: ClassVar[UnionTypeCPPSystem]
    _basis: BasisType
    _TypeBasis: type[BasisType]  # should by ClassVar, but cannot be nested yet

    def __init__(self, basis: BasisType) -> None:
        self._cpp = self._cpp_type(basis._cpp)  # type: ignore [reportPrivateUsage]
        self._update_basis()

    @classmethod
    def _from_cpp_object(cls: "type[Self]", cpp_obj: UnionCPPSystem) -> "Self":
        obj = cls.__new__(cls)
        obj._cpp = cpp_obj
        obj._update_basis()
        return obj

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.basis!r}, is_diagonal={self.is_diagonal})"

    def __str__(self) -> str:
        return self.__repr__().replace("Real", "").replace("Complex", "")

    def _update_basis(self) -> None:
        self._basis = self._TypeBasis._from_cpp_object(self._cpp.get_basis())  # type: ignore

    def diagonalize(
        self,
        diagonalizer: Diagonalizer = "eigen",
        float_type: FloatType = "float64",
        atol: float = 1e-6,
        sort_by_energy: bool = True,
        energy_range: tuple[Union["Quantity", None], Union["Quantity", None]] = (None, None),
        energy_unit: Optional[str] = None,
        m0: Optional[int] = None,
    ) -> "Self":
        """Diagonalize the Hamiltonian and update the basis to the eigenbasis.

        This method computes the eigenvalues and eigenvectors of the Hamiltonian
        and updates the internal basis of the system accordingly.

        Args:
            diagonalizer: The diagonalizer to use for the diagonalization.
                Possible values are "eigen", "lapacke_evd", "lapacke_evr", "feast",
                defaults to "eigen".
            float_type: The floating point precision to use for the diagonalization.
                Possible values are "float32", "float64", defaults to "float64".
            atol: The absolute tolerance, used to determine which entries of the sparse eigenbasis to keep.
                Defaults to 1e-6.
            sort_by_energy: Whether to sort the eigenstates by energy. Defaults to True.
            energy_range: A tuple specifying the energy range. Defaults to (None, None).
            energy_unit: The unit for the energy values. Defaults to None.
            m0: The search subspace size for the "feast" diagonalizer.

        Returns:
            Self: The updated instance of the system.

        """
        cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, self._cpp, float_type, m0=m0)

        min_energy_au, max_energy_au = energy_range
        if min_energy_au is not None:
            min_energy_au = QuantityScalar.from_pint_or_unit(min_energy_au, energy_unit, "ENERGY").to_base_unit()
        if max_energy_au is not None:
            max_energy_au = QuantityScalar.from_pint_or_unit(max_energy_au, energy_unit, "ENERGY").to_base_unit()
        self._cpp.diagonalize(cpp_diagonalizer, min_energy_au, max_energy_au, atol)

        if sort_by_energy:
            sorter = self._cpp.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            self._cpp.transform(sorter)

        self._update_basis()
        return self

    @property
    def is_diagonal(self) -> bool:
        return self._cpp.is_diagonal()

    @property
    def basis(self) -> BasisType:
        return self._basis

    @property
    def matrix(self) -> "csr_matrix":
        return self._cpp.get_matrix()

    def get_eigenbasis(self) -> BasisType:
        cpp_eigenbasis = self._cpp.get_eigenbasis()
        return self._TypeBasis._from_cpp_object(cpp_eigenbasis)  # type: ignore

    @overload
    def get_eigenvalues(self) -> "PlainQuantity[NDArray[Any]]": ...

    @overload
    def get_eigenvalues(self, unit: str) -> "NDArray[Any]": ...

    def get_eigenvalues(self, unit: Optional[str] = None):
        eigenvalues_au = self._cpp.get_eigenvalues()
        eigenvalues = QuantityArray.from_base_unit(eigenvalues_au, "ENERGY")
        return eigenvalues.to_pint_or_unit(unit)

    @overload
    def get_hamiltonian(self) -> "PlainQuantity[csr_matrix]": ...  # type: ignore [reportInvalidTypeArguments]

    @overload
    def get_hamiltonian(self, unit: str) -> "csr_matrix": ...

    def get_hamiltonian(self, unit: Optional[str] = None):
        hamiltonian_au = self._cpp.get_matrix()
        hamiltonian = QuantitySparse.from_base_unit(hamiltonian_au, "ENERGY")
        return hamiltonian.to_pint_or_unit(unit)


System = SystemBase[Any]
