# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from .utils import PairinteractionModule


@pytest.mark.parametrize("species", ["Yb171_mqdt", "Rb"])
def test_pair_potential(pi_module: PairinteractionModule, species: str) -> None:
    """Test multipole interaction."""
    ket = pi_module.KetAtom(species, nu=55.7, l=0, m=0.5)

    # Create a single-atom system
    basis = pi_module.BasisAtom(ket.species, nu=(ket.nu - 2, ket.nu + 2), l=(0, 3))
    print(f"Number of single-atom basis states: {basis.number_of_states}")

    system = pi_module.SystemAtom(basis)

    # Create two-atom systems for different interatomic distances and multipole orders
    delta_energy = 3  # GHz
    min_energy = 2 * ket.get_energy(unit="GHz") - delta_energy
    max_energy = 2 * ket.get_energy(unit="GHz") + delta_energy

    basis_pair = pi_module.BasisPair([system, system], energy=(min_energy, max_energy), energy_unit="GHz", m=(1, 1))
    print(f"Number of two-atom basis states: {basis_pair.number_of_states}")

    distances = np.linspace(0.2, 5, 5)
    system_pairs_0 = [pi_module.SystemPair(basis_pair) for d in distances]
    system_pairs_3 = [
        pi_module.SystemPair(basis_pair).set_interaction_order(3).set_distance(d, unit="micrometer") for d in distances
    ]
    system_pairs_4 = [
        pi_module.SystemPair(basis_pair).set_interaction_order(4).set_distance(d, unit="micrometer") for d in distances
    ]
    system_pairs_5 = [
        pi_module.SystemPair(basis_pair).set_interaction_order(5).set_distance(d, unit="micrometer") for d in distances
    ]

    # Separate the contributions of the different multipole orders
    order_3 = [
        a.get_hamiltonian(unit="GHz").todense() - b.get_hamiltonian(unit="GHz").todense()
        for a, b in zip(system_pairs_3, system_pairs_0, strict=True)
    ]
    order_4 = [
        a.get_hamiltonian(unit="GHz").todense() - b.get_hamiltonian(unit="GHz").todense()
        for a, b in zip(system_pairs_4, system_pairs_3, strict=True)
    ]
    order_5 = [
        a.get_hamiltonian(unit="GHz").todense() - b.get_hamiltonian(unit="GHz").todense()
        for a, b in zip(system_pairs_5, system_pairs_4, strict=True)
    ]

    # Check that each order of the multipole expansion of the interaction has a significant contribution
    # at short distance
    norm_3 = np.linalg.norm(order_3, axis=(1, 2))
    norm_4 = np.linalg.norm(order_4, axis=(1, 2))
    norm_5 = np.linalg.norm(order_5, axis=(1, 2))
    assert norm_3[0] * distances[0] ** 3 > 1
    assert norm_4[0] * distances[0] ** 4 > 1
    assert norm_5[0] * distances[0] ** 5 > 1

    # Check that for large/small distances, the lower/higher orders dominate
    assert norm_3[0] < norm_4[0] < norm_5[0]
    assert norm_3[-1] > norm_4[-1] > norm_5[-1]

    # Check that each order of the multipole expansion scales as expected
    assert np.allclose(norm_3 * distances**3, norm_3[0] * distances[0] ** 3)
    assert np.allclose(norm_4 * distances**4, norm_4[0] * distances[0] ** 4)
    assert np.allclose(norm_5 * distances**5, norm_5[0] * distances[0] ** 5)


def test_dipole_octupole_channel(pi_module: PairinteractionModule) -> None:
    """Test that pair states coupled only by dipole-octupole interaction appear at interaction order 5."""
    ket_s = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    ket_p = pi_module.KetAtom("Rb", n=60, l=1, j=1.5, m=1.5)
    ket_f = pi_module.KetAtom("Rb", n=58, l=3, j=3.5, m=-0.5)

    basis = pi_module.BasisAtom("Rb", additional_kets=[ket_s, ket_p, ket_f])
    system = pi_module.SystemAtom(basis)
    basis_pair = pi_module.BasisPair([system, system])

    # |s,s> and |p,f> differ by (delta_l1, delta_l2) = (1, 3), so they are coupled by the
    # dipole-octupole interaction but neither by dipole-dipole, dipole-quadrupole, nor
    # quadrupole-quadrupole interaction
    overlap_ss = basis_pair.get_overlaps([ket_s, ket_s])
    overlap_pf = basis_pair.get_overlaps([ket_p, ket_f])
    idx_ss = np.argmax(overlap_ss)
    idx_pf = np.argmax(overlap_pf)
    assert overlap_ss[idx_ss] > 0.99
    assert overlap_pf[idx_pf] > 0.99

    distances = [1, 2]  # micrometer
    elements = {}
    for order in [4, 5]:
        for distance in distances:
            system_pair = (
                pi_module.SystemPair(basis_pair).set_interaction_order(order).set_distance(distance, unit="micrometer")
            )
            hamiltonian = system_pair.get_hamiltonian(unit="GHz").toarray()
            elements[(order, distance)] = hamiltonian[idx_ss, idx_pf]

    # The coupling appears at order 5 and scales as 1/distance^5
    assert abs(elements[(4, 1)]) < 1e-10
    assert abs(elements[(5, 1)]) > 1e-5
    assert np.isclose(elements[(5, 1)], elements[(5, 2)] * 2**5)
