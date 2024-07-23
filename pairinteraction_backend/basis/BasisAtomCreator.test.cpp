#include "basis/BasisAtomCreator.hpp"

#include "basis/BasisAtom.hpp"
#include "database/Database.hpp"
#include "enums/Parity.hpp"
#include "enums/TransformationType.hpp"
#include "ket/KetAtom.hpp"
#include "ket/KetAtomCreator.hpp"
#include "utils/streamed.hpp"

#include <doctest/doctest.h>
#include <spdlog/spdlog.h>

DOCTEST_TEST_CASE("create a basis for strontium 88") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<float>()
                     .set_species("Sr88_singlet")
                     .restrict_quantum_number_n(60, 60)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    for (auto ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_singlet");
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*ket));
    }
}

DOCTEST_TEST_CASE("create a basis for strontium 87") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<float>()
                     .set_species("Sr87_mqdt")
                     .restrict_quantum_number_nu(59, 61)
                     .restrict_quantum_number_l(0, 0)
                     .create(database);
    for (auto ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr87_mqdt");
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*ket));
    }
}

DOCTEST_TEST_CASE("create a basis from kets") {
    Database &database = Database::get_global_instance();
    auto ket1 = KetAtomCreator<float>("Sr88_singlet", 59, 0, 0, 0).create(database);
    auto ket2 = KetAtomCreator<float>("Sr88_singlet", 60, 0, 0, 0).create(database);
    auto ket3 = KetAtomCreator<float>("Sr88_singlet", 61, 0, 0, 0).create(database);
    auto basis =
        BasisAtomCreator<float>().add_ket(ket1).add_ket(ket2).add_ket(ket3).create(database);
    for (auto ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_singlet");
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Ket: {}", fmt::streamed(*ket));
    }
}

DOCTEST_TEST_CASE("create a basis and sort it according to parity and m") {
    Database &database = Database::get_global_instance();
    auto basis_unsorted = BasisAtomCreator<float>()
                              .set_species("Rb")
                              .restrict_quantum_number_n(60, 60)
                              .restrict_quantum_number_l(0, 3)
                              .restrict_quantum_number_m(-0.5, 0.5)
                              .create(database);

    // Sort the basis by parity and the m quantum number
    auto basis = basis_unsorted->transformed(basis_unsorted->get_sorter(
        TransformationType::SORT_BY_PARITY | TransformationType::SORT_BY_QUANTUM_NUMBER_M));
    auto parity = Parity::ODD;
    auto quantum_number_m = std::numeric_limits<float>::lowest();
    for (size_t i = 0; i < basis->get_number_of_states(); ++i) {
        DOCTEST_CHECK(basis->get_parity(i) >= parity);
        if (basis->get_parity(i) != parity) {
            parity = basis->get_parity(i);
            quantum_number_m = std::numeric_limits<float>::lowest();
        }
        DOCTEST_CHECK(basis->get_quantum_number_m(i) >= quantum_number_m);
        quantum_number_m = basis->get_quantum_number_m(i);
    }

    // Check that the blocks are correctly determined
    auto blocks = basis->get_blocks(TransformationType::SORT_BY_PARITY |
                                    TransformationType::SORT_BY_QUANTUM_NUMBER_M);
    std::vector<size_t> expected_start = {0, 4, 8, 11};
    DOCTEST_CHECK(blocks.start.size() == expected_start.size());
    for (size_t i = 0; i < blocks.start.size(); ++i) {
        SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Block {} starts at {}", i, blocks.start[i]);
        DOCTEST_CHECK(blocks.start[i] == expected_start[i]);
    }

    // Test implicit conversion of an eigen matrix to a transformator
    size_t dim = basis->get_number_of_states();
    Eigen::SparseMatrix<float, Eigen::RowMajor> matrix(static_cast<long>(dim),
                                                       static_cast<long>(dim));
    matrix.setIdentity();
    auto transformed = basis->transformed(matrix);
    auto transformation = transformed->get_transformation();
    DOCTEST_CHECK(transformation.transformation_type.back() == TransformationType::ARBITRARY);
}
