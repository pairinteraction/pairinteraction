# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Any, Optional, Union, overload

from pairinteraction._wrapped.ket.ket_atom import KetAtom
from pairinteraction.perturbative.create_system import create_system_for_perturbative
from pairinteraction.perturbative.effective_hamiltonian import get_effective_hamiltonian_from_system
from pairinteraction.units import QuantityScalar

if TYPE_CHECKING:
    from pairinteraction._wrapped.ket.ket_pair import KetAtomTuple
    from pairinteraction._wrapped.system.system_pair import SystemPair
    from pairinteraction.units import PintArrayLike, PintFloat


@overload
def get_c6_from_system(ket_tuple: "KetAtomTuple", system_pair: "SystemPair[Any]", unit: None = None) -> "PintFloat": ...


@overload
def get_c6_from_system(ket_tuple: "KetAtomTuple", system_pair: "SystemPair[Any]", unit: str) -> float: ...


def get_c6_from_system(
    ket_tuple: "KetAtomTuple", system_pair: "SystemPair[Any]", unit: Optional[str] = None
) -> Union[float, "PintFloat"]:
    r"""Calculate the :math:`C_6` coefficient for a given pair state.

    This function takes a tuple of two KetAtom (i.e. (ketA, ketB) ),
    and calculates the :math:`C_6` coefficient of the corresponding pair state.
    If ketA and ketB are of the same species, they must be identical.
    If you want to calculate the second order perturbation corrections for two different states, of the same species,
    use the `get_effective_hamiltonian_from_system(...)` function instead.

    Args:
        ket_tuple: Tuple of two KetAtom: (ketA, ketB) for which the :math:`C_6` coefficient is calculated.
        system_pair: SystemPair object that defines the Hamiltonian
            (including electric and magnetic fields, the interatomic distance, ...)
        unit: The unit in which the :math:`C_6` coefficient will be returned.
            Default None will return a pint quantity.

    Returns:
        The :math:`C_6` coefficient.

    """
    if len(ket_tuple) != 2 or not isinstance(ket_tuple[0], KetAtom):
        raise ValueError("The C6 coefficient can only be calculated for a tuple of two KetAtoms as ket_tuple argument.")
    if ket_tuple[0].species == ket_tuple[1].species and ket_tuple[0] != ket_tuple[1]:
        raise ValueError(
            "If you want to calculate second order perturbation corrections for two different states "
            "use the `get_effective_hamiltonian_from_system(...)` function."
        )
    h_eff, _ = get_effective_hamiltonian_from_system(
        [ket_tuple], system_pair, perturbation_order=2, return_only_specified_order=True
    )
    c6_pint = h_eff[0, 0] * system_pair.get_distance() ** 6  # type: ignore [index] # PintArray does not know it can be indexed
    return QuantityScalar.from_pint(c6_pint, "c6").to_pint_or_unit(unit)


@overload
def get_c6(
    ket_tuple: "KetAtomTuple",
    distance_vector: Optional["PintArrayLike"] = None,
    electric_field: Optional["PintArrayLike"] = None,
    magnetic_field: Optional["PintArrayLike"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    *,
    unit: None = None,
) -> "PintFloat": ...


@overload
def get_c6(
    ket_tuple: "KetAtomTuple",
    distance_vector: Optional["PintArrayLike"] = None,
    electric_field: Optional["PintArrayLike"] = None,
    magnetic_field: Optional["PintArrayLike"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    *,
    unit: str,
) -> float: ...


def get_c6(
    ket_tuple: "KetAtomTuple",
    distance_vector: Optional["PintArrayLike"] = None,
    electric_field: Optional["PintArrayLike"] = None,
    magnetic_field: Optional["PintArrayLike"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    unit: Optional[str] = None,
) -> Union[float, "PintFloat"]:
    r"""Calculate the :math:`C_6` coefficient for a given tuple of ket states.

    This function calculates the :math:`C_6` coefficient in the desired unit. The input is a 2-tuple of single atom ket
    states.

    Args:
        ket_tuple: The input is a tuple repeating the same single atom state in the format (a,a).
        If a tuple with not exactly two identical states is given, a ValueError is raised.
        distance_vector: distance vector between the atoms.
        electric_field: electric field in the system.
        magnetic_field: magnetic field in the system.
        multipole_order: multipole-order of the interaction. Default is 3 (dipole-dipole).
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.
        unit: The unit to which to convert the result. Default None will return a pint quantity.

    Returns:
        The :math:`C_6` coefficient. If a unit is specified, the value in this unit is returned.

    """
    system_pair = create_system_for_perturbative(
        [ket_tuple],
        distance_vector,
        electric_field,
        magnetic_field,
        multipole_order,
        with_diamagnetism=with_diamagnetism,
        perturbation_order=2,
    )
    if distance_vector is None:
        system_pair.set_distance(100, unit="micrometer")
    return get_c6_from_system(ket_tuple, system_pair, unit=unit)
