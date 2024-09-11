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

    # TODO use "pint" and its contexts to handle units

    # Create systems for different values of the electric field
    systems = []
    for i in range(10):
        system = SystemAtomDouble(basis)
        system.set_electric_field([0, 0, i * 1e-4])
        systems.append(system)

    assert abs(systems[1].get_matrix().diagonal().sum() - systems[1].get_matrix().sum()) > 1e-10

    # TODO make system.set_electric_field, system.diagonalize, etc. return *this so that we can use
    # list comprehension to create the systems
    # TODO implement system.copy(), basis.copy(), etc.

    # Diagonalize the systems in parallel
    systems[1].diagonalize(diagonalizer, 12)

    assert abs(systems[1].get_matrix().diagonal().sum() - systems[1].get_matrix().sum()) < 1e-10

    # TODO support precision = 12 and no precision argument at all
    # TODO implement bindings to make diagonalize(systemms, diagonalizer) work
