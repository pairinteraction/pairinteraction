from collections.abc import Sequence
from typing import TYPE_CHECKING, Optional, TypeVar, Union

from pairinteraction import _backend
from pairinteraction._wrapped.cpp_types import (
    Diagonalizer,
    FloatType,
    get_cpp_diagonalize,
    get_cpp_diagonalizer,
)
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pint.facets.plain import PlainQuantity

    from pairinteraction._wrapped.system.System import System

    Quantity = TypeVar("Quantity", float, PlainQuantity[float])


def diagonalize(
    systems: Sequence["System"],
    diagonalizer: Diagonalizer = "eigen",
    float_type: FloatType = "float64",
    atol: float = 1e-6,
    sort_by_energy: bool = True,
    energy_range: tuple[Union["Quantity", None], Union["Quantity", None]] = (None, None),
    energy_unit: Optional[str] = None,
    m0: Optional[int] = None,
) -> None:
    """Diagonalize a list of systems in parallel using the C++ backend.

    A convenience function for diagonalizing a list of systems in parallel using the C++ backend.
    This is much faster than diagonalizing each system individually on the Python side.

    Examples:
        >>> import pairinteraction.real as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> systems = [pi.SystemAtom(basis).set_magnetic_field([0, 0, b], unit="gauss") for b in range(1, 4)]
        >>> print(systems[0])
        SystemAtom(BasisAtom(|Rb:58,S_1/2,-1/2⟩ ... |Rb:63,F_5/2,5/2⟩), is_diagonal=False)
        >>> pi.diagonalize(systems)
        >>> print(systems[0])
        SystemAtom(BasisAtom(|Rb:58,S_1/2,-1/2⟩ ... |Rb:63,F_5/2,5/2⟩), is_diagonal=True)

    Args:
        systems: A list of `SystemAtom` or `SystemPair` objects, which will get diagonalized inplace.
        diagonalizer: The diagonalizer method to use. Defaults to "eigen".
        float_type: The floating point precision to use for the diagonalization. Defaults to "float64".
        atol: The absolute tolerance allowed for eigenenergies and vectors
            Smaller values are set to 0 in the sparse eigenbasis. Defaults to 1e-6.
        sort_by_energy: Whether to sort the resulting basis by energy. Defaults to True.
        energy_range: A tuple specifying an energy range, in which eigenvlaues should be calculated.
            Specifying a range can speed up the diagonalization process (depending on the diagonalizer method).
            The accuracy of the eigenvalues is not affected by this, but not all eigenvalues will be calculated.
            Defaults to (None, None), i.e. calculate all eigenvalues.
        energy_unit: The unit in which the energy_range is given. Defaults to None assumes pint objects.
        m0: The search subspace size for the FEAST diagonalizer. Defaults to None.

    """
    cpp_systems = [s._cpp for s in systems]  # type: ignore [reportPrivateUsage]
    cpp_diagonalize_fct = get_cpp_diagonalize(systems[0])
    cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, cpp_systems[0], float_type, m0=m0)

    min_energy_au, max_energy_au = energy_range
    if min_energy_au is not None:
        min_energy_au = QuantityScalar.from_pint_or_unit(min_energy_au, energy_unit, "ENERGY").to_base_unit()
    if max_energy_au is not None:
        max_energy_au = QuantityScalar.from_pint_or_unit(max_energy_au, energy_unit, "ENERGY").to_base_unit()
    cpp_diagonalize_fct(cpp_systems, cpp_diagonalizer, min_energy_au, max_energy_au, atol)

    for system, cpp_system in zip(systems, cpp_systems):
        if sort_by_energy:
            sorter = cpp_system.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            cpp_system.transform(sorter)
        system._cpp = cpp_system  # type: ignore [reportPrivateUsage]
        system._update_basis()
