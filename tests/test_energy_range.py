"""Test energy range restriction in the diagonalization."""

import numpy as np

import pairinteraction.real as pi


def test_energy_range() -> None:
    """Test restricting the energy range in the diagonalization."""
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
    pair_energy = 2 * ket.get_energy(unit="GHz")
    distances = np.linspace(1, 4, 10)

    # Create a single-atom system
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 1))
    system = pi.SystemAtom(basis)

    # Create two-atom basis
    basis_pair = pi.BasisPair(
        [system, system], energy=(pair_energy - 10, pair_energy + 10), energy_unit="GHz", m=(1, 1)
    )

    # Diagonalize the systems for different distances in parallel and get all eigenenergies
    system_pairs = [pi.SystemPair(basis_pair).set_distance(d, unit="micrometer") for d in distances]
    pi.diagonalize(system_pairs, diagonalizer="eigen", sort_by_energy=True)
    eigenenergies_all = [system.get_eigenenergies(unit="GHz") for system in system_pairs]

    # Diagonalize the systems for different distances in parallel and get only the eigenenergies in an energy range
    system_pairs = [pi.SystemPair(basis_pair).set_distance(d, unit="micrometer") for d in distances]
    pi.diagonalize(
        system_pairs,
        diagonalizer="eigen",
        sort_by_energy=True,
        energy_range=(pair_energy - 5, pair_energy + 5),
        energy_unit="GHz",
    )
    eigenenergies_restricted = [system.get_eigenenergies(unit="GHz") for system in system_pairs]

    # Check the result
    eigenenergies_all_restricted = [
        eigenenergies[(eigenenergies < pair_energy + 5) & (eigenenergies > pair_energy - 5)]
        for eigenenergies in eigenenergies_all
    ]
    for e1, e2 in zip(eigenenergies_restricted, eigenenergies_all_restricted):
        np.testing.assert_allclose(e1, e2)
