// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetAtomCreator.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/ket/KetAtom.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("create a ket for rubidium") {
    Database &database = Database::get_global_instance();
    auto ket = KetAtomCreator("Rb", 60, 1, 0.5, 0.5).create(database);
    DOCTEST_CHECK(ket->get_species() == "Rb");
    DOCTEST_CHECK(ket->get_quantum_number("n") == 60);
    DOCTEST_CHECK(ket->get_quantum_number("l") == 1);
    DOCTEST_CHECK(ket->get_quantum_number("f") == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number("j") == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number("m") == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number("s") == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number("parity") == -1); // odd parity
    DOCTEST_MESSAGE("Ket: ", *ket);
}

DOCTEST_TEST_CASE("create a ket for strontium") {
    Database &database = Database::get_global_instance();
    auto ket = KetAtomCreator()
                   .set_species("Sr88_singlet")
                   .set_quantum_number("n", 60)
                   .set_quantum_number("l", 1)
                   .set_quantum_number("f", 1)
                   .set_quantum_number("m", 0)
                   .set_quantum_number("s", 0)
                   .create(database);
    DOCTEST_CHECK(ket->get_species() == "Sr88_singlet");
    DOCTEST_CHECK(ket->get_quantum_number("n") == 60);
    DOCTEST_CHECK(ket->get_quantum_number("f") == 1);
    DOCTEST_CHECK(ket->get_quantum_number("m") == 0);
    DOCTEST_CHECK(ket->get_quantum_number("parity") == -1); // odd parity
    DOCTEST_MESSAGE("Ket: ", *ket);
}

DOCTEST_TEST_CASE("quantum number std fallback") {
    Database &database = Database::get_global_instance();

    // n, f, and nu have no "std_<name>" column in the states table, so get_quantum_number_std
    // falls back to 0.
    auto ket = KetAtomCreator("Rb", 60, 0, 0.5, 0.5).create(database);
    DOCTEST_CHECK(ket->get_quantum_number_std("n") == 0);
    DOCTEST_CHECK(ket->get_quantum_number_std("f") == 0);
    DOCTEST_CHECK(ket->get_quantum_number_std("nu") == 0);

    // For an MQDT ket, quantum numbers that do have a "std_<name>" column report a nonzero spread,
    // confirming the zeros above come from the fallback and not from a stored value that is zero.
    auto ket_mqdt = KetAtomCreator()
                        .set_species("Yb171_mqdt")
                        .set_quantum_number("nu", 60)
                        .set_quantum_number("l", 0)
                        .set_quantum_number("f", 0.5)
                        .set_quantum_number("m", 0.5)
                        .create(database);
    DOCTEST_CHECK(ket_mqdt->get_quantum_number_std("l") > 0);
    DOCTEST_CHECK(ket_mqdt->get_quantum_number_std("n") == 0);
}

DOCTEST_TEST_CASE("test for equality") {
    Database &database = Database::get_global_instance();
    auto ket1 = KetAtomCreator("Rb", 60, 1, 0.5, 0.5).create(database);
    auto ket2 = KetAtomCreator("Rb", 60, 1, 0.5, 0.5).create(database);
    auto ket3 = KetAtomCreator("Rb", 60, 1, 1.5, 0.5).create(database);
    DOCTEST_CHECK(*ket1 == *ket2);
    DOCTEST_CHECK(*ket1 != *ket3);
}

} // namespace pairinteraction
