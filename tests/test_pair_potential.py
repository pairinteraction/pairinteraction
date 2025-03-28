"""Test the pair potential calculation."""

from pathlib import Path

import numpy as np
import pytest

import pairinteraction.real as pi

reference_kets_file = Path(__file__).parent.parent / "data/reference_pair_potential/kets.txt"
reference_eigenenergies_file = Path(__file__).parent.parent / "data/reference_pair_potential/eigenenergies.txt"
reference_overlaps_file = Path(__file__).parent.parent / "data/reference_pair_potential/overlaps.txt"


def test_pair_potential(generate_reference: bool) -> None:
    """Test calculating a pair potential."""
    # Create a single-atom system
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    print(f"Number of single-atom basis states: {basis.number_of_states}")

    system = pi.SystemAtom(basis)

    # Create two-atom systems for different interatomic distances
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
    delta_energy = 3  # GHz
    min_energy = 2 * ket.get_energy(unit="GHz") - delta_energy
    max_energy = 2 * ket.get_energy(unit="GHz") + delta_energy

    basis_pair = pi.BasisPair([system, system], energy=(min_energy, max_energy), energy_unit="GHz", m=(1, 1))
    print(f"Number of two-atom basis states: {basis_pair.number_of_states}")

    distances = np.linspace(1, 5, 5)
    system_pairs = [pi.SystemPair(basis_pair).set_distance(d, unit="micrometer") for d in distances]

    # Diagonalize the systems in parallel
    pi.diagonalize(system_pairs, diagonalizer="eigen", sort_by_energy=True)

    # Get the overlap with |ket, ket>
    overlaps = np.array([system.get_eigenbasis().get_overlaps([ket, ket]) for system in system_pairs])

    # Ensure that the overlaps sum up to one
    np.testing.assert_allclose(np.sum(overlaps, axis=1), np.ones(5))

    # Compare to reference data
    kets = [repr(ket) for ket in basis_pair.kets]
    eigenenergies = np.array([system.get_eigenenergies(unit="GHz") for system in system_pairs])
    eigenbasis = np.array([system.get_eigenbasis().get_coefficients().todense().A1 for system in system_pairs])

    if generate_reference:
        reference_kets_file.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_kets_file, kets, fmt="%s", delimiter="\t")
        np.savetxt(reference_eigenenergies_file, eigenenergies)
        np.savetxt(reference_overlaps_file, overlaps)
        pytest.skip("Reference data generated, skipping comparison test")

    np.testing.assert_equal(kets, np.loadtxt(reference_kets_file, dtype=str, delimiter="\t"))
    np.testing.assert_allclose(eigenenergies, np.loadtxt(reference_eigenenergies_file))
    np.testing.assert_allclose(overlaps, np.loadtxt(reference_overlaps_file), atol=1e-10)

    # Because of degeneracies, checking the eigenbasis against reference data is complicated.
    # Thus, we only check their normalization and orthogonality.
    cumulative_norm = (np.array(eigenbasis) * np.array(eigenbasis).conj()).sum(axis=1)
    np.testing.assert_allclose(cumulative_norm, 19 * np.ones(5))
