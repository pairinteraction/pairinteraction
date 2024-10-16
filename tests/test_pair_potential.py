"""
Test the Stark map calculation.
"""

from pathlib import Path

import numpy as np
import pytest
from pint import UnitRegistry

from pairinteraction.backend import (
    BasisAtomCreatorDouble,
    BasisCombinedCreatorDouble,
    Database,
    DiagonalizerEigenDouble,
    KetAtomCreatorDouble,
    SystemAtomDouble,
    SystemCombinedDouble,
    TransformationType,
    diagonalizeSystemCombinedDouble,
)

reference_kets_file = Path(__file__).parent.parent / "data/reference_pair_potential/kets.txt"
reference_eigenvalues_file = Path(__file__).parent.parent / "data/reference_pair_potential/eigenvalues.txt"


def test_pair_potential(
    ureg: UnitRegistry, generate_reference: bool, database_dir: str, download_missing: bool
) -> None:
    """Test calculating a pairpotential."""
    database = Database(download_missing, True, database_dir)
    diagonalizer = DiagonalizerEigenDouble()

    # Create a single-atom system
    basis = (
        BasisAtomCreatorDouble()
        .set_species("Rb")
        .restrict_quantum_number_n(58, 62)
        .restrict_quantum_number_l(0, 2)
        .create(database)
    )
    print(f"Number of single-atom basis states: {basis.get_number_of_states()}")

    system = SystemAtomDouble(basis)

    # Create two-atom systems for different interatomic distances
    ket = (
        KetAtomCreatorDouble()
        .set_species("Rb")
        .set_quantum_number_n(60)
        .set_quantum_number_l(0)
        .set_quantum_number_m(0.5)
        .create(database)
    )
    min_energy = 2 * ket.get_energy() * ureg.hartree - (3 * ureg.gigahertz).to(ureg.hartree, "spectroscopy")
    max_energy = 2 * ket.get_energy() * ureg.hartree + (3 * ureg.gigahertz).to(ureg.hartree, "spectroscopy")

    combined_basis = (
        BasisCombinedCreatorDouble()
        .add(system)
        .add(system)
        .restrict_energy(min_energy.to_base_units().magnitude, max_energy.to_base_units().magnitude)
        .restrict_quantum_number_m(1, 1)
        .create()
    )
    print(f"Number of two-atom basis states: {combined_basis.get_number_of_states()}")

    distances = np.linspace(1, 5, 5) * ureg.micrometer
    combined_systems = [
        SystemCombinedDouble(combined_basis).set_interatomic_distance(d.to_base_units().magnitude) for d in distances
    ]

    # Diagonalize the systems in parallel
    diagonalizeSystemCombinedDouble(combined_systems, diagonalizer)

    # Sort by the eigenvalues
    combined_systems = [
        system.transformed(system.get_sorter([TransformationType.SORT_BY_ENERGY])) for system in combined_systems
    ]

    # Compare to reference data
    kets = [str(ket) for ket in combined_basis.get_kets()]
    eigenvalues = [
        (system.get_eigenvalues() * ureg.hartree).to("gigahertz", "spectroscopy").magnitude
        for system in combined_systems
    ]
    eigenstates = [system.get_eigenstates().todense().A1 for system in combined_systems]

    if generate_reference:
        reference_kets_file.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_kets_file, kets, fmt="%s", delimiter=",")
        np.savetxt(reference_eigenvalues_file, eigenvalues)
        pytest.skip("Reference data generated, skipping comparison test")

    np.testing.assert_equal(kets, np.loadtxt(reference_kets_file, dtype=str, delimiter=","))
    np.testing.assert_allclose(eigenvalues, np.loadtxt(reference_eigenvalues_file))

    # Because of degeneracies, checking the eigenstates against reference data is complicated.
    # Thus, we only check their normalization and orthogonality.
    cumulative_norm = (np.array(eigenstates) * np.array(eigenstates).conj()).sum(axis=1)
    np.testing.assert_allclose(cumulative_norm, 19 * np.ones(5))
