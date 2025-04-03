# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test multipole interaction."""

import numpy as np
import pytest

import pairinteraction.real as pi


@pytest.mark.parametrize("species", ["Yb171_mqdt", "Rb"])
def test_pair_potential(species: str) -> None:
    """Test multipole interaction."""
    ket = pi.KetAtom(species, nu=55.7, l=0, m=0.5)

    # Create a single-atom system
    basis = pi.BasisAtom(ket.species, nu=(ket.nu - 2, ket.nu + 2), l=(0, 3))
    print(f"Number of single-atom basis states: {basis.number_of_states}")

    system = pi.SystemAtom(basis)

    # Create two-atom systems for different interatomic distances and multipole orders
    delta_energy = 3  # GHz
    min_energy = 2 * ket.get_energy(unit="GHz") - delta_energy
    max_energy = 2 * ket.get_energy(unit="GHz") + delta_energy

    basis_pair = pi.BasisPair([system, system], energy=(min_energy, max_energy), energy_unit="GHz", m=(1, 1))
    print(f"Number of two-atom basis states: {basis_pair.number_of_states}")

    distances = np.linspace(0.2, 5, 5)
    system_pairs_0 = [pi.SystemPair(basis_pair) for d in distances]
    system_pairs_3 = [pi.SystemPair(basis_pair).set_order(3).set_distance(d, unit="micrometer") for d in distances]
    system_pairs_4 = [pi.SystemPair(basis_pair).set_order(4).set_distance(d, unit="micrometer") for d in distances]

    # Separate the contributions of the different multipole orders
    order_3 = [
        a.get_hamiltonian(unit="GHz").todense() - b.get_hamiltonian(unit="GHz").todense()
        for a, b in zip(system_pairs_3, system_pairs_0)
    ]
    order_4 = [
        a.get_hamiltonian(unit="GHz").todense() - b.get_hamiltonian(unit="GHz").todense()
        for a, b in zip(system_pairs_4, system_pairs_3)
    ]

    # Check that each order of the multipole expansion of the interaction has a significant contribution
    # at short distance
    norm_3 = np.linalg.norm(order_3, axis=(1, 2))
    norm_4 = np.linalg.norm(order_4, axis=(1, 2))
    assert norm_3[0] * distances[0] ** 3 > 1
    assert norm_4[0] * distances[0] ** 4 > 1

    # Check that for large/small distances, dipole-dipole/dipole-quadrupole interaction dominates
    assert norm_3[0] < norm_4[0]
    assert norm_3[-1] > norm_4[-1]

    # Check that each order of the multipole expansion scales as expected
    assert np.allclose(norm_3 * distances**3, norm_3[0] * distances[0] ** 3)
    assert np.allclose(norm_4 * distances**4, norm_4[0] * distances[0] ** 4)
