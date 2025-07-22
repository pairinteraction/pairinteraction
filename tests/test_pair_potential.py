# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test the pair potential calculation."""

import numpy as np
import pairinteraction.real as pi
import pytest

from .compare_utils import REFERENCE_PATHS, compare_eigensystem_to_reference


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
    eigenvectors = np.array([system.get_eigenbasis().get_coefficients().todense().A1 for system in system_pairs])

    reference_path = REFERENCE_PATHS["pair_potential"]
    if generate_reference:
        reference_path.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_path / "kets.txt", kets, fmt="%s", delimiter="\t")
        np.savetxt(reference_path / "eigenenergies.txt", eigenenergies)
        np.savetxt(reference_path / "overlaps.txt", overlaps)
        pytest.skip("Reference data generated, skipping comparison test")

    compare_eigensystem_to_reference(reference_path, eigenenergies, overlaps, eigenvectors, kets)
