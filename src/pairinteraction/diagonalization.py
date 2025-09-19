# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later
from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Literal, TypeVar, Union, overload

from typing_extensions import deprecated

from pairinteraction import _backend
from pairinteraction.enums import get_cpp_float_type
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from pairinteraction.enums import FloatType
    from pairinteraction.system import SystemBase
    from pairinteraction.units import PintFloat

    Quantity = TypeVar("Quantity", bound=Union[float, "PintFloat"])


Diagonalizer = Literal["eigen", "lapacke_evd", "lapacke_evr", "feast"]
UnionCPPDiagonalizer = Union[_backend.DiagonalizerInterfaceReal, _backend.DiagonalizerInterfaceComplex]
UnionCPPDiagonalizerType = Union[type[_backend.DiagonalizerInterfaceReal], type[_backend.DiagonalizerInterfaceComplex]]

_DiagonalizerDict: dict[str, dict[Diagonalizer, UnionCPPDiagonalizerType]] = {
    "real": {
        "eigen": _backend.DiagonalizerEigenReal,
        "lapacke_evd": _backend.DiagonalizerLapackeEvdReal,
        "lapacke_evr": _backend.DiagonalizerLapackeEvrReal,
        "feast": _backend.DiagonalizerFeastReal,
    },
    "complex": {
        "eigen": _backend.DiagonalizerEigenComplex,
        "lapacke_evd": _backend.DiagonalizerLapackeEvdComplex,
        "lapacke_evr": _backend.DiagonalizerLapackeEvrComplex,
        "feast": _backend.DiagonalizerFeastComplex,
    },
}


@overload
def diagonalize(
    systems: Sequence[SystemBase[Any]],
    diagonalizer: Diagonalizer = "eigen",
    float_type: FloatType = "float64",
    rtol: float = 1e-6,
    sort_by_energy: bool = True,
    energy_range: tuple[Quantity | None, Quantity | None] = (None, None),
    energy_range_unit: str | None = None,
    m0: int | None = None,
) -> None: ...


@overload
@deprecated("Use energy_range_unit=... instead of energy_unit=...")
def diagonalize(
    systems: Sequence[SystemBase[Any]],
    diagonalizer: Diagonalizer = "eigen",
    float_type: FloatType = "float64",
    rtol: float = 1e-6,
    sort_by_energy: bool = True,
    energy_range: tuple[Quantity | None, Quantity | None] = (None, None),
    *,
    energy_unit: str | None,
    m0: int | None = None,
) -> None: ...


def diagonalize(
    systems: Sequence[SystemBase[Any]],
    diagonalizer: Diagonalizer = "eigen",
    float_type: FloatType = "float64",
    rtol: float = 1e-6,
    sort_by_energy: bool = True,
    energy_range: tuple[Quantity | None, Quantity | None] = (None, None),
    energy_range_unit: str | None = None,
    m0: int | None = None,
    *,
    energy_unit: str | None = None,
) -> None:
    """Diagonalize a list of systems in parallel using the C++ backend.

    A convenience function for diagonalizing a list of systems in parallel using the C++ backend.
    This is much faster than diagonalizing each system individually on the Python side.

    Examples:
        >>> import pairinteraction as pi
        >>> ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
        >>> basis = pi.BasisAtom("Rb", n=(58, 63), l=(0, 3))
        >>> systems = [pi.SystemAtom(basis).set_magnetic_field([0, 0, b], unit="gauss") for b in range(1, 4)]
        >>> print(systems[0])
        SystemAtom(BasisAtom(n=(58, 63), l=(0, 3)), is_diagonal=False)
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
        energy_range_unit: The unit in which the energy_range is given. Defaults to None assumes pint objects.
        m0: The search subspace size for the FEAST diagonalizer. Defaults to None.
        energy_unit: Deprecated, use energy_range_unit instead.

    """
    if energy_unit is not None:
        if energy_range_unit is not None:
            raise ValueError("energy_unit=... was replaced by energy_range_unit=..., do not use both together.")
        warnings.warn(
            "The energy_unit=... argument is deprecated, use energy_range_unit=... instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        energy_range_unit = energy_unit

    cpp_systems = [s._cpp for s in systems]
    cpp_diagonalize_fct = get_cpp_diagonalize_function(systems[0])
    cpp_diagonalizer = get_cpp_diagonalizer(diagonalizer, systems[0], float_type, m0=m0)

    energy_range_au: list[float | None] = [None, None]
    for i, energy in enumerate(energy_range):
        if energy is not None:
            energy_range_au[i] = QuantityScalar.convert_user_to_au(energy, energy_range_unit, "energy")
    cpp_diagonalize_fct(cpp_systems, cpp_diagonalizer, energy_range_au[0], energy_range_au[1], rtol)

    for system, cpp_system in zip(systems, cpp_systems):
        if sort_by_energy:
            sorter = cpp_system.get_sorter([_backend.TransformationType.SORT_BY_ENERGY])
            cpp_system.transform(sorter)
        system._cpp = cpp_system


def get_cpp_diagonalize_function(system: SystemBase[Any]) -> Callable[..., None]:
    if isinstance(system._cpp, _backend.SystemAtomReal):
        return _backend.diagonalizeSystemAtomReal
    if isinstance(system._cpp, _backend.SystemAtomComplex):
        return _backend.diagonalizeSystemAtomComplex
    if isinstance(system._cpp, _backend.SystemPairReal):
        return _backend.diagonalizeSystemPairReal
    if isinstance(system._cpp, _backend.SystemPairComplex):
        return _backend.diagonalizeSystemPairComplex
    raise TypeError(
        f"system must be of type SystemAtomReal, SystemPairReal, SystemAtomComplex, or SystemPairComplex, "
        f"not {type(system)}"
    )


def get_cpp_diagonalizer(
    diagonalizer: Diagonalizer,
    system: SystemBase[Any],
    float_type: FloatType,
    m0: int | None = None,
) -> UnionCPPDiagonalizer:
    if diagonalizer == "feast" and m0 is None:
        raise ValueError("m0 must be specified for the 'feast' diagonalizer")
    if diagonalizer != "feast" and m0 is not None:
        raise ValueError("m0 must not be specified if the diagonalizer is not 'feast'")

    if isinstance(system._cpp, (_backend.SystemAtomReal, _backend.SystemPairReal)):
        type_ = "real"
    elif isinstance(system._cpp, (_backend.SystemAtomComplex, _backend.SystemPairComplex)):
        type_ = "complex"
    else:
        raise TypeError(
            f"system must be of type SystemAtomReal, SystemPairReal, SystemAtomComplex, or SystemPairComplex, "
            f"not {type(system)}"
        )

    try:
        diagonalizer_class = _DiagonalizerDict[type_][diagonalizer]
    except KeyError:
        raise ValueError(
            f"Unknown diagonalizer '{diagonalizer}', should be one of {list(_DiagonalizerDict[type_].keys())}"
        ) from None

    cpp_float_type = get_cpp_float_type(float_type)
    if diagonalizer == "feast":
        return diagonalizer_class(m0=m0, float_type=cpp_float_type)  # type: ignore [call-arg]
    return diagonalizer_class(float_type=cpp_float_type)  # type: ignore [call-arg]
