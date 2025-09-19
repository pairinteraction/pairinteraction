# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING, Optional, Union, overload

import numpy as np

from pairinteraction.ket import KetAtom
from pairinteraction.state.state_base import StateBase

if TYPE_CHECKING:
    from typing_extensions import Self

    from pairinteraction.basis import BasisAtom
    from pairinteraction.database import Database
    from pairinteraction.enums import OperatorType
    from pairinteraction.units import PintComplex, PintFloat

logger = logging.getLogger(__name__)


class StateAtom(StateBase["BasisAtom", KetAtom]):
    """State of a single atom.

    A coefficient vector and a list of kets are used to represent an arbitrary single-atom state.

    Examples:
        >>> import pairinteraction as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(57, 63), l=(0, 3))
        >>> state = basis.get_corresponding_state(ket)
        >>> print(state)
        StateAtom(1.00 |Rb:60,S_1/2,1/2⟩)
        >>> ket2 = pi.KetAtom("Rb", n=60, l=1, j=0.5, m=0.5)
        >>> state2 = pi.StateAtom(ket2, basis)
        >>> print((2 * state2 - state).normalize())
        StateAtom(0.89 |Rb:60,P_1/2,1/2⟩ + -0.45 |Rb:60,S_1/2,1/2⟩)


    """

    _basis: "BasisAtom"
    _ket_class = KetAtom

    def __init__(self, ket: KetAtom, basis: "BasisAtom") -> None:
        """Initialize a state object representing a ket in a given basis.

        Args:
            ket: The ket to represent in the state.
            basis: The basis to which the state belongs.

        """
        new_basis: BasisAtom = basis.get_corresponding_state(ket)._basis
        ket_idx = new_basis.kets.index(ket)
        coeffs = new_basis._cpp.get_coefficients() * 0  # type: ignore [operator]
        coeffs[ket_idx, 0] = 1.0
        new_basis._cpp.set_coefficients(coeffs)
        self._basis = new_basis

    def __add__(self, other: "Self") -> "Self":
        """Add two states together.

        Args:
            other: The other state to add.

        Returns:
            A new state object representing the sum of the two states.

        """
        if not isinstance(other, type(self)):
            raise TypeError(f"Cannot add {type(self)} and {type(other)}.")
        if not all(ket_self == ket_other for ket_self, ket_other in zip(self.kets, other.kets)):
            raise ValueError("Cannot add states with different kets as basis.")

        coeffs = self._basis.get_coefficients() + other._basis.get_coefficients()
        new_basis = self._basis.copy()
        new_basis._cpp.set_coefficients(coeffs)
        return type(self)._from_basis_object(new_basis)

    def __sub__(self, other: "Self") -> "Self":
        """Subtract two states.

        Args:
            other: The other state to subtract.

        Returns:
            A new state object representing the difference of the two states.

        """
        return self.__add__(-1 * other)

    def __mul__(self, factor: complex) -> "Self":
        """Multiply the state with a scalar.

        Args:
            factor: The scalar to multiply with.

        Returns:
            A new state object representing the product of the state and the scalar.

        """
        if not isinstance(factor, (int, float, complex)):
            raise TypeError(f"Cannot multiply {type(self)} with {type(factor)}.")
        coeffs = factor * self._basis.get_coefficients()
        new_basis = self._basis.copy()
        new_basis._cpp.set_coefficients(coeffs)
        return type(self)._from_basis_object(new_basis)

    def __truediv__(self, factor: complex) -> "Self":
        """Divide the state by a scalar.

        Args:
            factor: The scalar to divide by.

        Returns:
            A new state object representing the quotient of the state and the scalar.

        """
        return self.__mul__(1 / factor)

    __rmul__ = __mul__  # for reverse multiplication, i.e. scalar * state will use state.__rmul__

    def normalize(self) -> "Self":
        """Normalize the coefficients of the state."""
        coeffs = self._basis.get_coefficients()
        self._basis._cpp.set_coefficients(coeffs / self.norm)
        return self

    def is_normalized(self, tol: float = 1e-10) -> bool:
        """Check if the state is normalized within a given tolerance.

        Args:
            tol: The tolerance for the normalization check. Default is 1e-10.

        Returns:
            True if the state is normalized within the given tolerance, False otherwise.

        """
        return abs(self.norm - 1) < tol  # type: ignore [return-value] # numpy

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
        return np.count_nonzero(self.get_coefficients()) == 1  # type: ignore [no-any-return]

    def get_amplitude(self, other: Union["Self", KetAtom]) -> Union[float, complex]:
        """Calculate the amplitude of the state with respect to another state or ket.

        This means the inner product <self|other>.

        Args:
            other: Either a state or a ket for which the amplitude should be calculated.

        Returns:
            The amplitude between self and other.

        """
        if not self.is_normalized() or (isinstance(other, StateAtom) and not other.is_normalized()):
            logger.warning("WARNING: get_amplitude is called with a non-normalized state.")
        return self._basis.get_amplitudes(other)[0]  # type: ignore [no-any-return]

    def get_overlap(self, other: Union["Self", KetAtom]) -> float:
        r"""Calculate the overlap of the state with respect to another state or ket.

        This means calculate :math:`|\langle \mathrm{self} | \mathrm{other} \rangle|^2`.

        Args:
            other: Either a state or a ket for which the overlap should be calculated.

        Returns:
            The overlap between self and other.

        """
        if not self.is_normalized() or (isinstance(other, StateAtom) and not other.is_normalized()):
            logger.warning("WARNING: get_overlap is called with a non-normalized state.")
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
        if not self.is_normalized() or (isinstance(other, StateAtom) and not other.is_normalized()):
            logger.warning("WARNING: get_matrix_element is called with a non-normalized state.")
        return self._basis.get_matrix_elements(other, operator, q, unit=unit)[0]  # type: ignore [index,no-any-return] # PintArray does not know it can be indexed


class StateAtomReal(StateAtom):
    _basis: "BasisAtom"
    _ket_class = KetAtom
