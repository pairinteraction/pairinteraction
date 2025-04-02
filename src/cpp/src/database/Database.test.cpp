// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

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

    AtomDescriptionByParameters description;
    description.quantum_number_n = 60;
    description.quantum_number_l = 0;
    description.quantum_number_m = 0.5;

    auto ket = database.get_ket("Rb", description);

    DOCTEST_MESSAGE("KetAtom: ", *ket);
}

DOCTEST_TEST_CASE("too large quantum number m") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByParameters description;
    description.quantum_number_n = 60;
    description.quantum_number_l = 0;
    description.quantum_number_m = 1.5;

    DOCTEST_CHECK_THROWS(database.get_ket("Rb", description));
}

DOCTEST_TEST_CASE("not uniquely specified ket") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByParameters description;
    description.quantum_number_n = 60;
    description.quantum_number_l = 0.9;
    description.quantum_number_m = 0.5;

    DOCTEST_CHECK_THROWS(database.get_ket("Rb", description));
}

DOCTEST_TEST_CASE("uniquely specified ket") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByParameters description;
    description.quantum_number_n = 60;
    description.quantum_number_l = 0.9;
    description.quantum_number_j = 0.5;
    description.quantum_number_m = 0.5;

    DOCTEST_CHECK_NOTHROW(database.get_ket("Rb", description));
}

DOCTEST_TEST_CASE("get a BasisAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges description;
    description.range_quantum_number_n = {60, 60};
    description.range_quantum_number_l = {0, 1};

    auto basis = database.get_basis<double>("Rb", description, {});

    for (const auto &ket : *basis) {
        DOCTEST_MESSAGE("KetAtom: ", *ket);
    }
}

DOCTEST_TEST_CASE("get an OperatorAtom") {
    Database &database = Database::get_global_instance();

    AtomDescriptionByRanges description;
    description.range_quantum_number_n = {60, 60};
    description.range_quantum_number_l = {0, 1};

    auto basis = database.get_basis<double>("Rb", description, {});

    auto dipole =
        database.get_matrix_elements<double>(basis, basis, OperatorType::ELECTRIC_DIPOLE, 0);

    DOCTEST_MESSAGE("Number of basis states: ", basis->get_number_of_states());
    DOCTEST_MESSAGE("Number of non-zero entries: ", dipole.nonZeros());
}
} // namespace pairinteraction
