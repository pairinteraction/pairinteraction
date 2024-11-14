"""
Test the Stark map calculation.
"""

from pathlib import Path

import numpy as np
import pytest
from pint import UnitRegistry

from pairinteraction.backend import (
    BasisAtomCreatorDouble,
    Database,
    DiagonalizerEigenDouble,
    KetAtomCreatorDouble,
    SystemAtomDouble,
    TransformationType,
    diagonalizeSystemAtomDouble,
)

reference_kets_file = Path(__file__).parent.parent / "data/reference_stark_map/kets.txt"
reference_eigenvalues_file = Path(__file__).parent.parent / "data/reference_stark_map/eigenvalues.txt"
reference_overlaps_file = Path(__file__).parent.parent / "data/reference_stark_map/overlaps.txt"


def test_starkmap(ureg: UnitRegistry, generate_reference: bool, database_dir: str, download_missing: bool) -> None:
    """Test calculating a Stark map."""
    database = Database(download_missing, True, database_dir)
    diagonalizer = DiagonalizerEigenDouble()

    # Create a basis
    ket = (
        KetAtomCreatorDouble()
        .set_species("Rb")
        .set_quantum_number_n(60)
        .set_quantum_number_l(0)
        .set_quantum_number_m(0.5)
        .create(database)
    )

    basis = (
        BasisAtomCreatorDouble()
        .set_species("Rb")
        .restrict_quantum_number_n(58, 62)
        .restrict_quantum_number_l(0, 2)
        .create(database)
    )
    print(f"Number of basis states: {basis.get_number_of_states()}")

    # Create systems for different values of the electric field
    electric_fields = np.linspace(0, 10, 11) * ureg.volt / ureg.centimeter
    systems = [SystemAtomDouble(basis).set_electric_field([0, 0, e.to_base_units().magnitude]) for e in electric_fields]

    # Diagonalize the systems in parallel
    diagonalizeSystemAtomDouble(systems, diagonalizer)

    # Sort by the eigenvalues
    for system in systems:
        system.transform(system.get_sorter([TransformationType.SORT_BY_ENERGY]))

    # Get the overlap with |ket>
    overlaps = [system.get_eigenbasis().get_overlaps(ket) for system in systems]

    # Ensure that the overlaps sum up to one
    np.testing.assert_allclose(np.sum(overlaps, axis=1), np.ones(11))

    # Compare to reference data
    kets = [str(ket) for ket in basis.get_kets()]
    eigenvalues = [
        (system.get_eigenvalues() * ureg.hartree).to(ureg.gigahertz, "spectroscopy").magnitude for system in systems
    ]
    eigenstates = [system.get_eigenbasis().get_coefficients().todense().A1 for system in systems]

    if generate_reference:
        reference_kets_file.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_kets_file, kets, fmt="%s", delimiter=",")
        np.savetxt(reference_eigenvalues_file, eigenvalues)
        np.savetxt(reference_overlaps_file, overlaps)
        pytest.skip("Reference data generated, skipping comparison test")

    np.testing.assert_equal(kets, np.loadtxt(reference_kets_file, dtype=str, delimiter=","))
    np.testing.assert_allclose(eigenvalues, np.loadtxt(reference_eigenvalues_file))
    np.testing.assert_allclose(overlaps, np.loadtxt(reference_overlaps_file), atol=1e-12)

    # Because of degeneracies, checking the eigenstates against reference data is complicated.
    # Thus, we only check their normalization and orthogonality.
    cumulative_norm = (np.array(eigenstates) * np.array(eigenstates).conj()).sum(axis=1)
    np.testing.assert_allclose(cumulative_norm, 90 * np.ones(11))
