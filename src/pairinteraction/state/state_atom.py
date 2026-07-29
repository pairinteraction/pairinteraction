# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import logging
import warnings
from typing import TYPE_CHECKING, cast, overload

import numpy as np
from scipy.sparse import csr_matrix
from typing_extensions import deprecated

from pairinteraction.enums import get_cpp_operator_type
from pairinteraction.ket import KetAtom, KetAtomReal
from pairinteraction.state.state_base import StateBase
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from collections.abc import Sequence

    from typing_extensions import Self

    from pairinteraction import _backend
    from pairinteraction.basis import BasisAtom
    from pairinteraction.database import Database
    from pairinteraction.enums import OperatorType
    from pairinteraction.units import PintComplex, PintFloat

logger = logging.getLogger(__name__)


class StateAtom(StateBase[KetAtom]):
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
        >>> state2 = ket2.to_state()
        >>> print((2 * state2 - state).normalize())
        StateAtom(0.89 |Rb:60,P_1/2,1/2⟩ - 0.45 |Rb:60,S_1/2,1/2⟩)
        >>> print(pi.StateAtom([2, 1], [ket, ket2]).normalize())
        StateAtom(0.89 |Rb:60,S_1/2,1/2⟩ + 0.45 |Rb:60,P_1/2,1/2⟩)


    """

    _cpp: _backend.BasisAtomComplex
    _ket_class = KetAtom

    @overload
    def __init__(
        self, coefficients: Sequence[complex], kets: Sequence[KetAtom], *, basis: BasisAtom | None = None
    ) -> None: ...

    @overload
    @deprecated("Use ket.to_state() instead of StateAtom(ket, basis).")
    def __init__(self, ket: KetAtom, basis: BasisAtom) -> None: ...

    def __init__(  # type: ignore [misc]
        self,
        coefficients: Sequence[complex] | KetAtom,
        kets: Sequence[KetAtom] | BasisAtom | None = None,
        *,
        basis: BasisAtom | None = None,
    ) -> None:
        """Initialize a state object from a coefficient vector and the corresponding kets.

        Args:
            coefficients: The coefficient of each of the given kets.
            kets: The kets the state is composed of.
            basis: The basis in which the state should be expressed.
                If None (default), a minimal basis consisting only of the given kets is constructed.
                Providing a basis is only relevant if you want the state to already live in a larger Hilbert space;
                when adding states, their bases are merged automatically.
                All given kets must be part of this basis.
                Since the coefficients are always defined with respect to the kets,
                only the kets of the given basis are used and the coefficients of the basis are ignored,
                i.e. the basis is canonicalized first.

        """
        super().__init__()

        if isinstance(coefficients, KetAtom):  # deprecated interface StateAtom(ket, basis)
            coefficients, kets, basis = self._unpack_deprecated_args(coefficients, kets, basis)
        if kets is None:
            raise TypeError("StateAtom.__init__() missing 1 required positional argument: 'kets'")
        kets = cast("Sequence[KetAtom]", kets)

        is_real = isinstance(self, StateAtomReal)
        coeffs = np.array(coefficients, dtype=float if is_real else complex).ravel()
        if len(coeffs) != len(kets):
            raise ValueError(f"Got {len(coeffs)} coefficients for {len(kets)} kets, these must match.")
        if len(kets) != len(set(kets)):
            raise ValueError("The given kets must be unique.")

        if basis is None:
            from pairinteraction.basis.basis_atom import get_cpp_basis_atom_from_kets

            cpp_basis = get_cpp_basis_atom_from_kets(kets, real=is_real)
        else:
            cpp_basis = basis._cpp.canonicalized()

        cpp_ket_to_index = {cpp_ket: i for i, cpp_ket in enumerate(cpp_basis.get_kets())}
        basis_coeffs = np.zeros((len(cpp_ket_to_index), 1), dtype=coeffs.dtype)
        for coeff, ket in zip(coeffs, kets, strict=True):
            if ket._cpp not in cpp_ket_to_index:
                raise ValueError(f"The ket {ket} is not part of the given basis.")
            basis_coeffs[cpp_ket_to_index[ket._cpp], 0] = coeff

        state_cpp = cpp_basis.get_state(0)  # single-state basis, i.e. shape (n_kets, 1)
        self._cpp = state_cpp.copy_with_coefficients(csr_matrix(basis_coeffs))

    @staticmethod
    def _unpack_deprecated_args(
        ket: KetAtom, kets: Sequence[KetAtom] | BasisAtom | None, basis: BasisAtom | None
    ) -> tuple[Sequence[complex], Sequence[KetAtom], BasisAtom]:
        """Translate the arguments of the deprecated StateAtom(ket, basis) interface into the new interface."""
        from pairinteraction.basis import BasisAtom  # imported here to avoid a circular import

        warnings.warn(
            "Calling StateAtom(ket, basis) is deprecated use ket.to_state() instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        if kets is not None:  # in the deprecated interface the second argument was the basis
            if basis is not None:
                raise TypeError("The basis must not be given both as positional and as keyword argument.")
            basis = cast("BasisAtom", kets)
        if not isinstance(basis, BasisAtom):
            raise TypeError("The basis must be given as a BasisAtom object when creating a StateAtom from a ket.")
        return [1], [ket], basis

    def __add__(self, other: Self | KetAtom) -> Self:
        """Add two states together.

        The bases are merged into a common basis containing the kets of both states,
        and the coefficients are re-expressed in this merged basis before adding.

        Args:
            other: The other state to add. A ket is converted to a state via `ket.to_state()`.

        Returns:
            A new state object representing the sum of the two states.

        """
        if isinstance(other, KetAtom):
            other = cast("Self", other.to_state())
        if type(self) is not type(other):
            raise TypeError(f"Cannot add/subtract {type(self)} and {type(other)}.")

        # merge the (canonical) bases and re-express the coefficients in the merged basis;
        merged_cpp = self._cpp.canonicalized().merge(other._cpp.canonicalized())
        cpp_op = get_cpp_operator_type("identity")
        coeffs1 = self._cpp.get_matrix_elements(merged_cpp, cpp_op, 0)
        coeffs2 = other._cpp.get_matrix_elements(merged_cpp, cpp_op, 0)
        coeffs = coeffs1 + coeffs2
        new_cpp = merged_cpp.get_state(0)  # single-state basis, i.e. shape (n_kets, 1)

        new_cpp = new_cpp.copy_with_coefficients(coeffs)
        return type(self)._from_cpp_object(new_cpp)

    def __sub__(self, other: Self | KetAtom) -> Self:
        """Subtract two states.

        Args:
            other: The other state to subtract. A ket is converted to a state via `ket.to_state()`.

        Returns:
            A new state object representing the difference of the two states.

        """
        if isinstance(other, KetAtom):
            other = cast("Self", other.to_state())
        return self.__add__(-1 * other)

    def __mul__(self, factor: complex) -> Self:
        """Multiply the state with a scalar.

        Args:
            factor: The scalar to multiply with.

        Returns:
            A new state object representing the product of the state and the scalar.

        """
        if not isinstance(factor, (int, float, complex)):
            raise TypeError(f"Cannot multiply {type(self)} with {type(factor)}.")
        coeffs = factor * self._cpp.get_coefficients()  # type: ignore [operator]
        new_cpp = self._cpp.copy_with_coefficients(coeffs)
        return type(self)._from_cpp_object(new_cpp)

    def __truediv__(self, factor: complex) -> Self:
        """Divide the state by a scalar.

        Args:
            factor: The scalar to divide by.

        Returns:
            A new state object representing the quotient of the state and the scalar.

        """
        return self.__mul__(1 / factor)

    def __neg__(self) -> Self:
        """Negate the state.

        Returns:
            A new state object with all coefficients multiplied by minus one.

        """
        return self.__mul__(-1)

    __rmul__ = __mul__  # for reverse multiplication, i.e. scalar * state will use state.__rmul__

    def normalize(self) -> Self:
        """Normalize the coefficients of the state."""
        coeffs = self._cpp.get_coefficients()
        self._cpp = self._cpp.copy_with_coefficients(coeffs / self.norm)
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
    def database(self) -> Database:
        """The database used for this object."""
        return self.get_ket(0).database

    @property
    def species(self) -> str:
        """The atomic species."""
        return self.get_ket(0).species

    @property
    def is_canonical(self) -> bool:
        return np.count_nonzero(self.get_coefficients()) == 1  # type: ignore [no-any-return]

    def get_amplitude(self, other: Self | KetAtom) -> float | complex:
        """Calculate the amplitude of the state with respect to another state or ket.

        This means the inner product <self|other>.

        Args:
            other: Either a state or a ket for which the amplitude should be calculated.

        Returns:
            The amplitude between self and other.

        """
        return self.get_matrix_element(other, "identity", 0, unit="")

    def get_overlap(self, other: Self | KetAtom) -> float:
        r"""Calculate the overlap of the state with respect to another state or ket.

        This means calculate :math:`|\langle \mathrm{self} | \mathrm{other} \rangle|^2`.

        Args:
            other: Either a state or a ket for which the overlap should be calculated.

        Returns:
            The overlap between self and other.

        """
        return abs(self.get_amplitude(other)) ** 2

    @overload
    def get_matrix_element(
        self, other: KetAtom | Self, operator: OperatorType, q: int, unit: None = None
    ) -> PintFloat | PintComplex: ...

    @overload
    def get_matrix_element(
        self, other: KetAtom | Self, operator: OperatorType, q: int, unit: str
    ) -> float | complex: ...

    def get_matrix_element(
        self, other: KetAtom | Self, operator: OperatorType, q: int, unit: str | None = None
    ) -> PintFloat | PintComplex | float | complex:
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
            logger.warning("get_matrix_element/get_overlap/get_amplitude is called with a non-normalized state.")

        cpp_op = get_cpp_operator_type(operator)

        if isinstance(other, KetAtom):
            other = cast("Self", other.to_state())
        if isinstance(other, StateAtom):
            matrix_elements_au = self._cpp.get_matrix_elements(other._cpp, cpp_op, q).toarray().ravel()[0]
            return QuantityScalar.convert_au_to_user(matrix_elements_au, operator, unit)
        raise TypeError(f"Unknown type: {type(other)=}")


class StateAtomReal(StateAtom):
    _cpp: _backend.BasisAtomReal  # type: ignore [assignment]
    _ket_class = KetAtomReal
