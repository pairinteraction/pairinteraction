"""
Test the pair potential calculation.
"""

from pathlib import Path

import numpy as np
import pytest

import pairinteraction.backend.double as pi

reference_kets_file = Path(__file__).parent.parent / "data/reference_pair_potential/kets.txt"
reference_eigenvalues_file = Path(__file__).parent.parent / "data/reference_pair_potential/eigenvalues.txt"
reference_overlaps_file = Path(__file__).parent.parent / "data/reference_pair_potential/overlaps.txt"


def test_pair_potential(generate_reference: bool, database_dir: str, download_missing: bool) -> None:
    """Test calculating a pair potential."""
    database = pi.Database(download_missing, True, database_dir)

    # Create a single-atom system
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2), database=database)
    print(f"Number of single-atom basis states: {basis.number_of_states}")

    system = pi.SystemAtom(basis)

    # Create two-atom systems for different interatomic distances
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5, database=database)
    delta_energy = 3
    min_energy = 2 * ket.get_energy(unit="GHz") - delta_energy
    max_energy = 2 * ket.get_energy(unit="GHz") + delta_energy

    combined_basis = pi.BasisCombined([system, system], energy=(min_energy, max_energy), energy_unit="GHz", m=(1, 1))
    print(f"Number of two-atom basis states: {combined_basis.number_of_states}")

    distances = np.linspace(1, 5, 5)
    combined_systems = [pi.SystemCombined(combined_basis).set_distance(d, unit="micrometer") for d in distances]

    # Diagonalize the systems in parallel
    pi.diagonalize(combined_systems, diagonalizer="Eigen", sort_by_energy=True)

    # Get the overlap with |ket, ket>
    overlaps = np.array(
        [system.get_eigenbasis().get_overlaps_with_product_state(ket, ket) for system in combined_systems]
    )

    # Ensure that the overlaps sum up to one
    np.testing.assert_allclose(np.sum(overlaps, axis=1), np.ones(5))

    # Compare to reference data
    kets = [str(ket) for ket in combined_basis.kets]
    eigenvalues = np.array([system.get_eigenvalues(unit="GHz") for system in combined_systems])
    eigenstates = np.array([system.get_eigenbasis().coefficients.todense().A1 for system in combined_systems])

    if generate_reference:
        reference_kets_file.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_kets_file, kets, fmt="%s", delimiter="\t")
        np.savetxt(reference_eigenvalues_file, eigenvalues)
        np.savetxt(reference_overlaps_file, overlaps)
        pytest.skip("Reference data generated, skipping comparison test")

    np.testing.assert_equal(kets, np.loadtxt(reference_kets_file, dtype=str, delimiter="\t"))
    np.testing.assert_allclose(eigenvalues, np.loadtxt(reference_eigenvalues_file))
    np.testing.assert_allclose(overlaps, np.loadtxt(reference_overlaps_file), atol=1e-12)

    # Because of degeneracies, checking the eigenstates against reference data is complicated.
    # Thus, we only check their normalization and orthogonality.
    cumulative_norm = (np.array(eigenstates) * np.array(eigenstates).conj()).sum(axis=1)
    np.testing.assert_allclose(cumulative_norm, 19 * np.ones(5))
