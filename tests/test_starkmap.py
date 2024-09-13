from pathlib import Path

from pairinteraction.backend import BasisAtomCreatorDouble, Database, DiagonalizerEigenDouble, SystemAtomDouble

# TODO from pairinteraction.backend_double import BasisAtomCreator, Database, DiagonalizerEigen, SystemAtom

databasedir = Path(__file__).parent.parent / "data/database"


def test_starkmap():
    database = Database(False, True, databasedir)
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

    # Create systems for different values of the electric field
    systems = [SystemAtomDouble(basis).set_electric_field([0, 0, i * 1e-4]) for i in range(10)]

    assert abs(systems[1].get_matrix().diagonal().sum() - systems[1].get_matrix().sum()) > 1e-10

    # Diagonalize the systems in parallel
    systems[1].diagonalize(diagonalizer)

    # TODO diagonalize(systemms, diagonalizer)

    assert abs(systems[1].get_matrix().diagonal().sum() - systems[1].get_matrix().sum()) < 1e-10
