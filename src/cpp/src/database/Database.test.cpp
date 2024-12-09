#include "pairinteraction/database/Database.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/AtomDescriptionByParameters.hpp"
#include "pairinteraction/database/AtomDescriptionByRanges.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("get a KetAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByParameters<float> description;
    description.quantum_number_n = 60;
    description.quantum_number_l = 0;
    description.quantum_number_m = 0;

    auto ket = database.get_ket<float>("Rb", description);

    DOCTEST_MESSAGE("KetAtom: ", *ket);
}

DOCTEST_TEST_CASE("get a BasisAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges<float> description;
    description.range_quantum_number_n = {60, 60};
    description.range_quantum_number_l = {0, 1};

    auto basis = database.get_basis<float>("Rb", description, {});

    for (const auto &ket : *basis) {
        DOCTEST_MESSAGE("KetAtom: ", *ket);
    }
}

DOCTEST_TEST_CASE("get an OperatorAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges<float> description;
    description.range_quantum_number_n = {60, 60};
    description.range_quantum_number_l = {0, 1};

    auto basis = database.get_basis<float>("Rb", description, {});

    auto dipole =
        database.get_matrix_elements<float>(basis, basis, OperatorType::ELECTRIC_DIPOLE, 0);

    DOCTEST_MESSAGE("Number of basis states: ", basis->get_number_of_states());
    DOCTEST_MESSAGE("Number of non-zero entries: ", dipole.nonZeros());
}
} // namespace pairinteraction
