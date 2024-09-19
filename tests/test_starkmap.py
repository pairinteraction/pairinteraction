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

reference_data_file = Path(__file__).parent.parent / "data/reference_stark_map/rb_small_basis.txt"


def test_starkmap(ureg: UnitRegistry, generate_reference: bool, database_dir: Path, download_missing: bool) -> None:
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

    electric_fields = np.linspace(0, 10, 10) * ureg.volt / ureg.centimeter

    # Create systems for different values of the electric field
    systems = [SystemAtomDouble(basis).set_electric_field([0, 0, e.to_base_units().magnitude]) for e in electric_fields]

    # Diagonalize the systems in parallel
    diagonalize(systems, diagonalizer)
    eigenvalues = [
        (system.get_matrix().diagonal().real * ureg.hartree).to("gigahertz", "spectroscopy") for system in systems
    ]

    # Compare eigenvalues to reference values
    if generate_reference:
        reference_data_file.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(reference_data_file, [calc.magnitude for calc in eigenvalues])
        print(f"Reference data generated and saved to {reference_data_file}")
        pytest.skip("Reference data generated, skipping comparison test")

    reference_values = np.loadtxt(reference_data_file)
    for calc, ref in zip(eigenvalues, reference_values):
        np.testing.assert_allclose(calc.magnitude, ref)
