"""
Test energy range restriction in the diagonalization.
"""

import numpy as np

import pairinteraction.backend.double as pi


def test_energy_range(database_dir: str, download_missing: bool) -> None:
    """Test restricting the energy range in the diagonalization."""
    pi.initialize_global_database(download_missing, True, database_dir)

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

    # Diagonalize the systems for different distances in parallel and get all eigenvalues
    system_pairs = [pi.SystemPair(basis_pair).set_distance(d, unit="micrometer") for d in distances]
    pi.diagonalize(system_pairs, diagonalizer="Eigen", sort_by_energy=True)
    eigenvalues_all = [system.get_eigenvalues(unit="GHz") for system in system_pairs]

    # Diagonalize the systems for different distances in parallel and get only the eigenvalues in an energy range
    system_pairs = [pi.SystemPair(basis_pair).set_distance(d, unit="micrometer") for d in distances]
    pi.diagonalize(
        system_pairs,
        diagonalizer="Eigen",
        sort_by_energy=True,
        energy_range=(pair_energy - 5, pair_energy + 5),
        energy_unit="GHz",
    )
    eigenvalues_restricted = [system.get_eigenvalues(unit="GHz") for system in system_pairs]

    # Check the result
    eigenvalues_all_restricted = [
        eigenvalues[(eigenvalues < pair_energy + 5) & (eigenvalues > pair_energy - 5)]
        for eigenvalues in eigenvalues_all
    ]
    for e1, e2 in zip(eigenvalues_restricted, eigenvalues_all_restricted):
        np.testing.assert_allclose(e1, e2)
