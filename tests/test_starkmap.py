# SPDX-FileCopyrightText: 2024 PairInteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from .utils import REFERENCE_PATHS, compare_eigensystem_to_reference

if TYPE_CHECKING:
    from .utils import PairinteractionModule


def test_starkmap(pi_module: PairinteractionModule, generate_reference: bool) -> None:
    """Test calculating a Stark map."""
    # Create a basis
    ket = pi_module.KetAtom("Rb", n=60, l=0, m=0.5)
    basis = pi_module.BasisAtom("Rb", n=(58, 62), l=(0, 2))
    print(f"Number of basis states: {basis.number_of_states}")

    electric_fields = np.linspace(0, 10, 11)
    # Create systems for different values of the electric field
    systems = [pi_module.SystemAtom(basis).set_electric_field([0, 0, e], unit="V/cm") for e in electric_fields]

    # Diagonalize the systems in parallel
    pi_module.diagonalize(systems, diagonalizer="eigen", sort_by_energy=True)

    # Get the overlap with |ket>
    overlaps = np.array([system.basis.get_overlaps(ket) for system in systems])

    # Compare to reference data
    kets = [ket.get_label("raw") for ket in systems[0].basis.kets]
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
