# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from collections.abc import Collection
from typing import TYPE_CHECKING, Any, Optional, Union

import numpy as np

from pairinteraction import (
    complex as pi_complex,
    real as pi_real,
)
from pairinteraction.units import BaseUnits, QuantityScalar, ureg

if TYPE_CHECKING:
    from pairinteraction._wrapped.ket.ket_pair import KetAtomTuple
    from pairinteraction._wrapped.system.system_pair import SystemPair
    from pairinteraction.units import PintArrayLike, PintFloat

logger = logging.getLogger(__name__)


def create_system_for_perturbative(  # noqa: C901, PLR0912, PLR0915
    ket_tuples: Collection["KetAtomTuple"],
    distance_vector: Optional["PintArrayLike"] = None,
    electric_field: Optional["PintArrayLike"] = None,
    magnetic_field: Optional["PintArrayLike"] = None,
    multipole_order: int = 3,
    with_diamagnetism: bool = False,
    perturbation_order: int = 2,
    min_population_admixture: float = 1e-4,
) -> "SystemPair[Any]":
    r"""Create a good estimate for a system in which to perform perturbative calculations.

    This function takes a list of 2-tuples of ket states and creates a pair system holding a larger basis.
    The parameters of the basis are adjusted by the electric and magnetic field vectors of the system, as well
    as the distance and the multipole-order of the interaction. For higher-order perturbation theory, larger
    systems can be considered. Diamagnetism can be considered as well.

    Args:
        ket_tuples: List of all pair states that span up the model space. The system is created such that
        the effective Hamiltonian of the model system can be calculated accurately at a later stage.
        distance_vector: distance vector between the atoms.
        electric_field: electric field in the system.
        magnetic_field: magnetic field in the system.
        multipole_order: multipole-order of the interaction. Default is 3 (dipole-dipole).
        with_diamagnetism: True if diamagnetic term should be considered. Default is False.
        perturbation_order: order of perturbative calculation the system shall be used for. Default is 2.
        min_population_admixture: Minimum population admixture to still consider. This is used to estimate the
            delta energy that needs to be considered. Default is 1e-4.

    Returns:
        Pair system that can be used for perturbative calculations.

    """
    electric_field = electric_field if electric_field is not None else [0, 0, 0]
    magnetic_field = magnetic_field if magnetic_field is not None else [0, 0, 0]

    pi = pi_real if electric_field[1] == 0 and magnetic_field[1] == 0 else pi_complex  # type: ignore [index]
    is_cylindrical_symmetric = all(x == 0 for x in [*magnetic_field[:2], *electric_field[:2]])  # type: ignore [index]

    system_atoms: list[Union[pi_real.SystemAtom, pi_complex.SystemAtom]] = []
    max_dipoles_au: list[PintFloat] = []

    delta_n = 7
    delta_l = perturbation_order * (multipole_order - 2)
    for i in range(2):
        kets = [ket_tuple[i] for ket_tuple in ket_tuples]
        nlfm = np.transpose([[ket.n, ket.l, ket.f, ket.m] for ket in kets])
        n_range = (int(np.min(nlfm[0])) - delta_n, int(np.max(nlfm[0])) + delta_n)
        l_range = None
        if not any(ket.is_calculated_with_mqdt for ket in kets):
            l_range = (np.min(nlfm[1]) - delta_l, np.max(nlfm[1]) + delta_l)
        f_range = None
        m_range = None
        if is_cylindrical_symmetric:
            f_range = (np.min(nlfm[2]) - delta_l, np.max(nlfm[2]) + delta_l)
            m_range = (np.min(nlfm[3]) - delta_l, np.max(nlfm[3]) + delta_l)
        basis = pi.BasisAtom(kets[0].species, n=n_range, l=l_range, f=f_range, m=m_range)
        system = pi.SystemAtom(basis)
        system.set_diamagnetism_enabled(with_diamagnetism)
        system.set_magnetic_field(magnetic_field)
        system.set_electric_field(electric_field)
        system_atoms.append(system)

        electric_dipoles_au = [
            basis.get_matrix_elements(ket, "electric_dipole", q=q).to_base_units().magnitude
            for ket in kets
            for q in [+1, -1, 0]
        ]
        max_dipoles_au.append(np.max(np.abs(electric_dipoles_au)))

    pi.diagonalize(system_atoms)

    # if no distance vector is given, we set delta_energy_au to very small and use the min_number_of_kets below
    delta_energy_au = QuantityScalar.from_unit(0.1, "MHz", "energy").to_base_unit()

    if distance_vector is not None:
        # estimate max interaction energy from dipole-dipole interaction at given distance
        distance_vector_au = [0 if d == 0 else d.to_base_units().magnitude for d in distance_vector]  # type: ignore [union-attr]
        distance_au = np.linalg.norm(distance_vector_au)
        interaction_energy_max_au = (
            max_dipoles_au[0]
            * max_dipoles_au[1]
            * ureg.Quantity(1, ureg.coulomb_constant).to_base_units().magnitude
            / distance_au**3
        )

        # estimate delta energy we need to consider, to take into account perturbations from states
        # that are at least populated (according to first order perturbation theory) by min_population_admixture
        delta_energy_au = max(interaction_energy_max_au / np.sqrt(min_population_admixture), delta_energy_au)

    pair_energies_au = [
        sum(
            system.get_corresponding_energy(ket).to_base_units().magnitude
            for system, ket in zip(system_atoms, ket_tuple)
        )
        for ket_tuple in ket_tuples
    ]

    def get_basis_pair(delta_energy_au: float) -> Union[pi_real.BasisPair, pi_complex.BasisPair]:
        return pi.BasisPair(
            system_atoms,
            energy=(min(pair_energies_au) - delta_energy_au, max(pair_energies_au) + delta_energy_au),
            energy_unit=str(BaseUnits["energy"]),
        )

    basis_pair = get_basis_pair(delta_energy_au)

    # minimum number of kets in the pair basis, even if the estimate of delta energy would lead to a smaller basis
    # this is to avoid misleading results when considering large distances
    min_number_of_kets = 2_000
    if basis_pair.number_of_kets < min_number_of_kets and perturbation_order > 0:
        logger.debug(
            "The basis of the pair system estimated from the interaction energy is very small (%d kets)."
            " This might lead to misleading results, thus we will increase delta energy,"
            " to at least include %d number of states.",
            basis_pair.number_of_kets,
            min_number_of_kets,
        )
        min_delta, max_delta = delta_energy_au, None
        delta_energy_au = QuantityScalar.from_unit(1, "GHz", "energy").to_base_unit()

        # make a bisect search to get a sensible basis size between: min_number_of_kets and 1.5 * min_number_of_kets
        while delta_energy_au < 1:  # stop if delta_energy_au is 1 and the basis is still very small
            basis_pair = get_basis_pair(delta_energy_au)
            if basis_pair.number_of_kets < min_number_of_kets:
                min_delta = delta_energy_au
                if max_delta is None:
                    delta_energy_au *= 2
                else:
                    delta_energy_au = (delta_energy_au + max_delta) / 2
            elif basis_pair.number_of_kets > min_number_of_kets * 1.5:
                max_delta = delta_energy_au
                delta_energy_au = (delta_energy_au + min_delta) / 2
            else:
                break

    logger.debug("The pair basis for the perturbative calculations consists of %d kets.", basis_pair.number_of_kets)
    if basis_pair.number_of_kets > 20_000:
        logger.warning(
            "The basis of the pair system is very large (%d kets). "
            "This might lead to long calculation times. "
            "Consider using a smaller basis by adjusting the min_population_admixture.",
            basis_pair.number_of_kets,
        )

    system_pair = pi.SystemPair(basis_pair)
    if distance_vector is not None:
        system_pair.set_distance_vector(distance_vector)
    system_pair.set_interaction_order(multipole_order)
    return system_pair
