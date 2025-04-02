# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Any, Optional, TypeVar, Union

from pairinteraction import _backend, _wrapped
from pairinteraction._wrapped.diagonalize.diagonalizer import Diagonalizer, get_cpp_diagonalizer
from pairinteraction._wrapped.enums import FloatType
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pairinteraction._wrapped.system.system import SystemBase
    from pairinteraction.units import PintFloat

    Quantity = TypeVar("Quantity", bound=Union[float, "PintFloat"])


def diagonalize(
    systems: Sequence["SystemBase[Any]"],
    diagonalizer: Diagonalizer = "eigen",
    float_type: FloatType = "float64",
    rtol: float = 1e-6,
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
        rtol: The relative tolerance allowed for eigenenergies. The error in eigenenergies is bounded
            by rtol * ||H||, where ||H|| is the norm of the Hamiltonian matrix. Defaults to 1e-6.
        sort_by_energy: Whether to sort the resulting basis by energy. Defaults to True.
        energy_range: A tuple specifying an energy range, in which eigenvlaues should be calculated.
            Specifying a range can speed up the diagonalization process (depending on the diagonalizer method).
            The accuracy of the eigenenergies is not affected by this, but not all eigenenergies will be calculated.
            Defaults to (None, None), i.e. calculate all eigenenergies.
        energy_unit: The unit in which the energy_range is given. Defaults to None assumes pint objects.
        m0: The search subspace size for the FEAST diagonalizer. Defaults to None.

    """
    cpp_systems = [s._cpp for s in systems]
    cpp_diagonalize_fct = get_cpp_diagonalize_function(systems[0])
    cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, systems[0], float_type, m0=m0)

    energy_range_au: list[Optional[float]] = [None, None]
    for i, energy in enumerate(energy_range):
        if energy is not None:
            energy_range_au[i] = QuantityScalar.from_pint_or_unit(energy, energy_unit, "energy").to_base_unit()
    cpp_diagonalize_fct(cpp_systems, cpp_diagonalizer, energy_range_au[0], energy_range_au[1], rtol)

    for system, cpp_system in zip(systems, cpp_systems):
        if sort_by_energy:
            sorter = cpp_system.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            cpp_system.transform(sorter)
        system._cpp = cpp_system
        system._update_basis()


def get_cpp_diagonalize_function(system: "SystemBase[Any]") -> Callable[..., None]:
    if isinstance(system, _wrapped.SystemAtomReal):
        return _backend.diagonalizeSystemAtomReal
    if isinstance(system, _wrapped.SystemAtomComplex):
        return _backend.diagonalizeSystemAtomComplex
    if isinstance(system, _wrapped.SystemPairReal):
        return _backend.diagonalizeSystemPairReal
    if isinstance(system, _wrapped.SystemPairComplex):
        return _backend.diagonalizeSystemPairComplex
    raise TypeError(
        f"system must be of type SystemAtomReal, SystemPairReal, SystemAtomComplex, or SystemPairComplex, "
        f"not {type(system)}"
    )
