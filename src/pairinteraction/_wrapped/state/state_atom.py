# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, Optional, TypeVar, Union, overload

import numpy as np

from pairinteraction._wrapped.ket.ket_atom import KetAtom
from pairinteraction._wrapped.state.state import StateBase

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction._wrapped.basis.basis_atom import BasisAtom, BasisAtomComplex, BasisAtomReal
    from pairinteraction._wrapped.database.database import Database
    from pairinteraction._wrapped.enums import OperatorType
    from pairinteraction.units import PintComplex, PintFloat

BasisType = TypeVar("BasisType", bound="BasisAtom[Any]", covariant=True)


class StateAtom(StateBase[BasisType, KetAtom]):
    """State of a single atom.

    A coefficient vector and a list of kets are used to represent an arbitrary single atom state.

    Examples:
        >>> import pairinteraction.real as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(57, 63), l=(0, 3))
        >>> state = basis.get_corresponding_state(ket)
        >>> print(state)
        StateAtom(1.00 |Rb:60,S_1/2,1/2⟩)

    """

    @property
    def database(self) -> "Database":
        """The database used for this object."""
        return self._basis.database

    @property
    def species(self) -> str:
        """The atomic species."""
        return self.kets[0].species

    @property
    def is_canonical(self) -> bool:
        return np.count_nonzero(self.get_coefficients()) == 1

    def get_amplitude(self, other: Union["Self", KetAtom]) -> Union[float, complex]:
        """Calculate the amplitude of the state with respect to another state or ket.

        This means the inner product <self|other>.

        Args:
            other: Either a state or a ket for which the amplitude should be calculated.

        Returns:
            The amplitude between self and other.

        """
        return self._basis.get_amplitudes(other)[0]  # type: ignore [no-any-return]

    def get_overlap(self, other: Union["Self", KetAtom]) -> Union[float, complex]:
        """Calculate the overlap of the state with respect to another state or ket.

        This means the inner product |<self|other>|^2.

        Args:
            other: Either a state or a ket for which the overlap should be calculated.

        Returns:
            The overlap between self and other.

        """
        return self._basis.get_overlaps(other)[0]  # type: ignore [no-any-return]

    @overload
    def get_matrix_element(
        self, other: Union["KetAtom", "Self"], operator: "OperatorType", q: int, unit: None = None
    ) -> Union["PintFloat", "PintComplex"]: ...  # type: ignore [type-var] # see "PintComplex"

    @overload
    def get_matrix_element(
        self, other: Union["KetAtom", "Self"], operator: "OperatorType", q: int, unit: str
    ) -> Union[float, complex]: ...

    def get_matrix_element(
        self, other: Union["KetAtom", "Self"], operator: "OperatorType", q: int, unit: Optional[str] = None
    ) -> Union["PintFloat", "PintComplex", float, complex]:
        """Calculate the matrix element of the operator with respect to the state and another state or ket.

        This means the inner product <self|operator|other>.

        Args:
            other: Either a state or a ket for which the matrix element should be calculated.
            operator: The operator for which the matrix element should be calculated.
            q: The projection quantum number of the operator.
            unit: The unit in which the result should be returned.
                Default None will return a `pint.Quantity`.

        Returns:
            The matrix element between self and other.

        """
        return self._basis.get_matrix_elements(other, operator, q, unit=unit)[0]  # type: ignore [index,no-any-return] # PintArray does not know it can be indexed


class StateAtomReal(StateAtom["BasisAtomReal"]):
    _basis: "BasisAtomReal"
    _TypeKet = KetAtom


class StateAtomComplex(StateAtom["BasisAtomComplex"]):
    _basis: "BasisAtomComplex"
    _TypeKet = KetAtom
