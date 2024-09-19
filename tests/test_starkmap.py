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
    SystemAtomDouble,
    diagonalize,
)

reference_kets_file = Path(__file__).parent.parent / "data/reference_stark_map/kets.txt"
reference_eigenvalues_file = Path(__file__).parent.parent / "data/reference_stark_map/eigenvalues.txt"
# reference_eigenstates_file = Path(__file__).parent.parent / "data/reference_stark_map/eigenstates.txt"


def test_starkmap(ureg: UnitRegistry, generate_reference: bool, database_dir: str, download_missing: bool) -> None:
    """Test calculating a Stark map."""
    database = Database(download_missing, True, database_dir)
    diagonalizer = DiagonalizerEigenDouble()

    # Create a basis
    basis = (
        BasisAtomCreatorDouble()
        .set_species("Rb")
        .restrict_quantum_number_n(58, 62)
        .restrict_quantum_number_l(0, 2)
        .restrict_quantum_number_m(0.5, 0.5)
        .create(database)
    )
    print(f"Number of basis states: {basis.get_number_of_states()}")

    electric_fields = np.linspace(0, 10, 11) * ureg.volt / ureg.centimeter

    # Create systems for different values of the electric field
    systems = [SystemAtomDouble(basis).set_electric_field([0, 0, e.to_base_units().magnitude]) for e in electric_fields]

    # Diagonalize the systems in parallel
    diagonalize(systems, diagonalizer)

    # Compare to reference data
    kets = [str(ket) for ket in systems[0].get_basis().get_kets()]
    eigenvalues = [
        (system.get_matrix().diagonal().real * ureg.hartree).to("gigahertz", "spectroscopy").magnitude
        for system in systems
    ]
    # eigenstates = [system.get_basis().get_coefficients().todense().A1 for system in systems]

    if generate_reference:
        reference_kets_file.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_kets_file, kets, fmt="%s", delimiter=";")
        np.savetxt(reference_eigenvalues_file, eigenvalues)
        # np.savetxt(reference_eigenstates_file, eigenstates)
        pytest.skip("Reference data generated, skipping comparison test")

    np.testing.assert_equal(kets, np.loadtxt(reference_kets_file, dtype=str, delimiter=";"))
    np.testing.assert_allclose(eigenvalues, np.loadtxt(reference_eigenvalues_file))

    # Because of degeneracies, checking the eigenstates will require a more sophisticated comparison
    # than the simple approach below that is currently commented out.
    # np.testing.assert_allclose(eigenstates, np.loadtxt(reference_eigenstates_file))
