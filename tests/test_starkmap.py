# SPDX-FileCopyrightText: 2024 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

"""Test the Stark map calculation."""

import numpy as np
import pairinteraction.real as pi
import pytest

from .compare_utils import REFERENCE_PATHS, compare_eigensystem_to_reference


def test_starkmap(generate_reference: bool) -> None:
    """Test calculating a Stark map."""
    # Create a basis
    ket = pi.KetAtom("Rb", n=60, l=0, m=0.5)
    basis = pi.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    print(f"Number of basis states: {basis.number_of_states}")

    electric_fields = np.linspace(0, 10, 11)
    # Create systems for different values of the electric field
    systems = [pi.SystemAtom(basis).set_electric_field([0, 0, e], unit="V/cm") for e in electric_fields]

    # Diagonalize the systems in parallel
    pi.diagonalize(systems, diagonalizer="eigen", sort_by_energy=True)

    # Get the overlap with |ket>
    overlaps = np.array([system.basis.get_overlaps(ket) for system in systems])

    # Compare to reference data
    kets = [repr(ket) for ket in systems[0].basis.kets]
    eigenenergies = np.array([system.get_eigenenergies(unit="GHz") for system in systems])
    eigenvectors = np.array([system.get_eigenbasis().get_coefficients().todense().A1 for system in systems])

    reference_path = REFERENCE_PATHS["stark_map"]
    if generate_reference:
        reference_path.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_path / "kets.txt", kets, fmt="%s", delimiter="\t")
        np.savetxt(reference_path / "eigenenergies.txt", eigenenergies)
        np.savetxt(reference_path / "overlaps.txt", overlaps)
        pytest.skip("Reference data generated, skipping comparison test")

    compare_eigensystem_to_reference(reference_path, eigenenergies, overlaps, eigenvectors, kets)
