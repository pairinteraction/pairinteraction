#include "pintr/ket/KetAtomCreator.hpp"

#include "pintr/database/Database.hpp"
#include "pintr/ket/KetAtom.hpp"
#include "pintr/utils/streamed.hpp"

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

DOCTEST_TEST_CASE("create a ket for rubidium") {
    Database &database = Database::get_global_instance();
    auto ket = KetAtomCreator<float>("Rb", 60, 1, 0.5, 0.5).create(database);
    DOCTEST_CHECK(ket->get_species() == "Rb");
    DOCTEST_CHECK(ket->get_quantum_number_n() == 60);
    DOCTEST_CHECK(ket->get_quantum_number_l() == 1);
    DOCTEST_CHECK(ket->get_quantum_number_f() == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number_j() == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number_m() == 0.5);
    DOCTEST_CHECK(ket->get_quantum_number_s() == 0.5);
    DOCTEST_CHECK(ket->get_parity() == Parity::ODD);
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*ket));
}

DOCTEST_TEST_CASE("create a ket for strontium") {
    Database &database = Database::get_global_instance();
    auto ket = KetAtomCreator<float>()
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
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*ket));
}
