// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
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
    DOCTEST_CHECK(ket->get_quantum_number_n() == 60);
    DOCTEST_CHECK(ket->get_quantum_number_l() == 1);
    DOCTEST_CHECK(ket->get_quantum_number_f() == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number_j() == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number_m() == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number_s() == 0.5);
    DOCTEST_CHECK(ket->get_parity() == Parity::ODD);
    DOCTEST_MESSAGE("Ket: ", *ket);
}

DOCTEST_TEST_CASE("create a ket for strontium") {
    Database &database = Database::get_global_instance();
    auto ket = KetAtomCreator()
                   .set_species("Sr88_singlet")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(1)
                   .set_quantum_number_f(1)
                   .set_quantum_number_m(0)
                   .set_quantum_number_s(0)
                   .create(database);
    DOCTEST_CHECK(ket->get_species() == "Sr88_singlet");
    DOCTEST_CHECK(ket->get_quantum_number_n() == 60);
    DOCTEST_CHECK(ket->get_quantum_number_f() == 1);
    DOCTEST_CHECK(ket->get_quantum_number_m() == 0);
    DOCTEST_CHECK(ket->get_parity() == Parity::ODD);
    DOCTEST_MESSAGE("Ket: ", *ket);
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
