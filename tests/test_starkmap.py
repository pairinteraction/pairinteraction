"""
Test the Stark map calculation.
"""

from pathlib import Path

import numpy as np
import pytest

import pairinteraction.backend.double as pi

reference_kets_file = Path(__file__).parent.parent / "data/reference_stark_map/kets.txt"
reference_eigenvalues_file = Path(__file__).parent.parent / "data/reference_stark_map/eigenvalues.txt"
reference_overlaps_file = Path(__file__).parent.parent / "data/reference_stark_map/overlaps.txt"


def test_starkmap(generate_reference: bool, database_dir: str, download_missing: bool) -> None:
    """Test calculating a Stark map."""
    database = pi.Database(download_missing, True, database_dir)

    # Create a basis
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5, database=database)
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2), database=database)
    print(f"Number of basis states: {basis.number_of_states}")

    electric_fields = np.linspace(0, 10, 11)
    # Create systems for different values of the electric field
    systems = [pi.SystemAtom(basis).set_electric_field([0, 0, e], unit="V/cm") for e in electric_fields]

    # Diagonalize the systems in parallel
    pi.diagonalize(systems, diagonalizer="Eigen", sort_by_energy=True)

    # Get the overlap with |ket>
    overlaps = np.array([system.basis.get_overlaps(ket) for system in systems])

    # Ensure that the overlaps sum up to one
    np.testing.assert_allclose(np.sum(overlaps, axis=1), np.ones(len(electric_fields)))

    # Compare to reference data
    kets = [str(ket) for ket in systems[0].basis.kets]
    eigenvalues = np.array([system.get_eigenvalues(unit="GHz") for system in systems])
    eigenstates = np.array([system.get_eigenbasis().coefficients.todense().A1 for system in systems])

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
    np.testing.assert_allclose(cumulative_norm, 90 * np.ones(len(electric_fields)))
