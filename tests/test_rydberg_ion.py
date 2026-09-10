# SPDX-FileCopyrightText: 2025 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Tests for Rydberg ions, i.e., charged species whose interaction includes monopole terms."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from pairinteraction import _backend
from pairinteraction.units import ureg

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_monopole_matrix_elements(pi_module: PairinteractionModule) -> None:
    """The monopole operator is the total charge in units of the charge -e of the Rydberg electron."""
    ket_ion = pi_module.KetAtom("Sr88_ion", n=60, l=0, j=0.5, m=0.5)
    assert ket_ion.get_matrix_element(ket_ion, "electric_monopole", q=0, unit="e") == pytest.approx(-1)

    # The monopole operator is diagonal and the same for all states of the ion
    basis_ion = pi_module.BasisAtom("Sr88_ion", n=(59, 61), l=(0, 2))
    monopole = basis_ion.get_matrix_elements(basis_ion, "electric_monopole", q=0, unit="e").toarray()
    np.testing.assert_allclose(monopole, -np.eye(basis_ion.number_of_states))


@pytest.mark.parametrize("order", [2, 3])
@pytest.mark.parametrize("direction", ["z", "x", "y"])
def test_atom_ion_pair_vs_point_charge(
    pi_module: PairinteractionModule, use_real: bool, order: int, direction: str
) -> None:
    """A Rydberg ion in a single state must act on a Rydberg atom like a classical point charge."""
    if direction == "y" and use_real:
        pytest.skip("a y-component requires complex numbers")

    distance = 3  # micrometer
    distance_vector = {"x": [distance, 0, 0], "y": [0, distance, 0], "z": [0, 0, distance]}[direction]

    basis_atom = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 3), m=(0.5, 0.5) if direction == "z" else None)
    ket_ion = pi_module.KetAtom("Sr88_ion", n=60, l=0, j=0.5, m=0.5)
    basis_ion = pi_module.BasisAtom("Sr88_ion", n=(0, 0), additional_kets=[ket_ion])
    assert basis_ion.number_of_states == 1

    # Pair system of the atom and the ion
    basis_pair = pi_module.BasisPair([pi_module.SystemAtom(basis_atom), pi_module.SystemAtom(basis_ion)])
    assert basis_pair.number_of_states == basis_atom.number_of_states
    system_pair = (
        pi_module.SystemPair(basis_pair)
        .set_interaction_order(order)
        .set_distance_vector(distance_vector, unit="micrometer")
        .diagonalize(diagonalizer="eigen")
    )
    energies_pair = system_pair.get_eigenenergies(unit="GHz") - ket_ion.get_energy(unit="GHz")

    # Atom in the field of a classical point charge
    system_reference = (
        pi_module.SystemAtom(basis_atom)
        .set_ion_charge(1, unit="e")
        .set_ion_interaction_order(order)
        .set_ion_distance_vector(distance_vector, unit="micrometer")
        .diagonalize(diagonalizer="eigen")
    )
    energies_reference = system_reference.get_eigenenergies(unit="GHz")

    np.testing.assert_allclose(energies_pair, energies_reference, atol=1e-6, rtol=0)

    # The interaction must have a significant effect
    energies_unperturbed = pi_module.SystemAtom(basis_atom).diagonalize(diagonalizer="eigen").get_eigenenergies("GHz")
    assert np.linalg.norm(energies_pair - energies_unperturbed) > 1e-3


def test_ion_ion_coulomb_repulsion(pi_module: PairinteractionModule) -> None:
    """Two ions in a single state each interact via the repulsive Coulomb interaction Z1*Z2/R."""
    ket_ion = pi_module.KetAtom("Sr88_ion", n=60, l=0, j=0.5, m=0.5)
    basis_ion = pi_module.BasisAtom("Sr88_ion", n=(0, 0), additional_kets=[ket_ion])
    system_ion = pi_module.SystemAtom(basis_ion)
    basis_pair = pi_module.BasisPair([system_ion, system_ion])
    assert basis_pair.number_of_states == 1

    distances = np.array([0.5, 1, 2, 5])  # micrometer
    for order in [1, 2, 3]:
        energies = np.array(
            [
                pi_module.SystemPair(basis_pair)
                .set_interaction_order(order)
                .set_distance(d, unit="micrometer")
                .get_hamiltonian(unit="hartree")[0, 0]
                for d in distances
            ]
        )
        energies -= 2 * ket_ion.get_energy(unit="hartree")
        distances_au = ureg.Quantity(distances, "micrometer").to("bohr").magnitude
        np.testing.assert_allclose(energies, 1 / distances_au, rtol=1e-10)
        assert np.all(energies > 0)


def test_ion_ion_multipole_orders(pi_module: PairinteractionModule) -> None:
    """The contributions of the individual multipole orders scale as 1/R^order."""
    ket = pi_module.KetAtom("Sr88_ion", n=60, l=0, j=0.5, m=0.5)
    basis = pi_module.BasisAtom("Sr88_ion", n=(ket.n - 1, ket.n + 1), l=(0, 2))
    system = pi_module.SystemAtom(basis)

    delta_energy = 300  # GHz
    pair_energy = 2 * ket.get_energy(unit="GHz")
    basis_pair = pi_module.BasisPair(
        [system, system], energy=(pair_energy - delta_energy, pair_energy + delta_energy), energy_unit="GHz", m=(1, 1)
    )
    assert basis_pair.number_of_states > 1

    distances = np.linspace(1, 5, 5)
    hamiltonian_0 = pi_module.SystemPair(basis_pair).get_hamiltonian(unit="GHz").toarray()
    hamiltonians = {0: [hamiltonian_0 for _ in distances]}
    for order in [1, 2, 3]:
        hamiltonians[order] = [
            pi_module.SystemPair(basis_pair)
            .set_interaction_order(order)
            .set_distance(d, unit="micrometer")
            .get_hamiltonian(unit="GHz")
            .toarray()
            for d in distances
        ]
    contributions = {
        order: np.linalg.norm(np.array(hamiltonians[order]) - np.array(hamiltonians[order - 1]), axis=(1, 2))
        for order in [1, 2, 3]
    }

    for order, norm in contributions.items():
        assert norm[0] > 0
        np.testing.assert_allclose(norm * distances**order, norm[0] * distances[0] ** order, rtol=1e-8)

    # The monopole-monopole interaction is a constant shift of all pair states
    shift = hamiltonians[1][0] - hamiltonians[0][0]
    np.testing.assert_allclose(shift, shift[0, 0] * np.eye(basis_pair.number_of_states), atol=1e-10)


def test_neutral_pair_unaffected(pi_module: PairinteractionModule) -> None:
    """The monopole terms vanish for neutral atoms and must not alter the block structure of the Hamiltonian."""
    basis = pi_module.BasisAtom("Rb", n=(59, 61), l=(0, 2))
    system = pi_module.SystemAtom(basis)
    ket = pi_module.KetAtom("Rb", n=60, l=0, j=0.5, m=0.5)
    pair_energy = 2 * ket.get_energy(unit="GHz")
    basis_pair = pi_module.BasisPair([system, system], energy=(pair_energy - 3, pair_energy + 3), energy_unit="GHz")

    # For a distance vector along z, the Hamiltonian must be block-diagonal with respect to the quantum number m
    system_pair = pi_module.SystemPair(basis_pair).set_interaction_order(3).set_distance(3, unit="micrometer")
    labels = [_backend.SorterType.QUANTUM_NUMBER_M]
    system_pair._cpp.transform(system_pair._cpp.get_sorter(labels))
    blocks = system_pair._cpp.get_indices_of_blocks(labels)
    assert len(blocks) > 1
    hamiltonian = system_pair.get_hamiltonian(unit="GHz").toarray()
    for block in blocks:
        hamiltonian[block.start : block.end, block.start : block.end] = 0
    assert np.all(hamiltonian == 0)

    # Without charged species, the interaction orders 1 and 2 do not contribute
    system_pair_order_2 = pi_module.SystemPair(basis_pair).set_interaction_order(2).set_distance(3, unit="micrometer")
    system_pair_order_0 = pi_module.SystemPair(basis_pair)
    hamiltonian_order_2 = system_pair_order_2.get_hamiltonian(unit="GHz").toarray()
    hamiltonian_order_0 = system_pair_order_0.get_hamiltonian(unit="GHz").toarray()
    np.testing.assert_allclose(hamiltonian_order_2, hamiltonian_order_0)
