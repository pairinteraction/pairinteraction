#include "pairinteraction/basis/BasisAtomCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("create a basis for strontium 88") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<float>()
                     .set_species("Sr88_singlet")
                     .restrict_quantum_number_n(60, 60)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_singlet");
        DOCTEST_MESSAGE("Ket: ", *ket);
    }
}

DOCTEST_TEST_CASE("create a basis for strontium 87") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<float>()
                     .set_species("Sr87_mqdt")
                     .restrict_quantum_number_nu(59, 61)
                     .restrict_quantum_number_l(0, 0)
                     .create(database);
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr87_mqdt");
        DOCTEST_MESSAGE("Ket: ", *ket);
    }
}

DOCTEST_TEST_CASE("create a basis from kets") {
    Database &database = Database::get_global_instance();
    auto ket1 = KetAtomCreator<float>("Sr88_singlet", 59, 0, 0, 0).create(database);
    auto ket2 = KetAtomCreator<float>("Sr88_singlet", 60, 0, 0, 0).create(database);
    auto ket3 = KetAtomCreator<float>("Sr88_singlet", 61, 0, 0, 0).create(database);
    auto basis =
        BasisAtomCreator<float>().append_ket(ket1).append_ket(ket2).append_ket(ket3).create(
            database);
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_singlet");
        DOCTEST_MESSAGE("Ket: ", *ket);
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
    auto sorter = basis_unsorted->get_sorter(
        {TransformationType::SORT_BY_PARITY, TransformationType::SORT_BY_QUANTUM_NUMBER_M});
    auto basis = basis_unsorted->transformed(sorter);

    // Check if the basis is properly sorted
    auto parity = Parity::ODD;
    auto quantum_number_m = std::numeric_limits<float>::lowest();
    for (size_t i = 0; i < basis->get_number_of_states(); ++i) {
        DOCTEST_MESSAGE("State ", i, ": Parity = ", basis->get_parity(i),
                        ", M = ", basis->get_quantum_number_m(i));
        DOCTEST_CHECK(basis->get_parity(i) >= parity);
        if (basis->get_parity(i) != parity) {
            parity = basis->get_parity(i);
            quantum_number_m = std::numeric_limits<float>::lowest();
        }
        DOCTEST_CHECK(basis->get_quantum_number_m(i) >= quantum_number_m);
        quantum_number_m = basis->get_quantum_number_m(i);
    }

    // Check that the blocks are correctly determined
    auto blocks = basis->get_indices_of_blocks(
        {TransformationType::SORT_BY_PARITY, TransformationType::SORT_BY_QUANTUM_NUMBER_M});
    std::vector<size_t> expected_start = {0, 4, 8, 11};

    DOCTEST_CHECK(blocks.size() == expected_start.size());

    size_t idx = 0;
    for (const auto &block : blocks) {
        DOCTEST_MESSAGE("Block ", idx, " starts at ", block.start);
        DOCTEST_CHECK(block.start == expected_start[idx]);
        idx++;
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

DOCTEST_TEST_CASE("calculation of matrix elements") {
    auto &database = Database::get_global_instance();

    auto ket_s = KetAtomCreator<double>()
                     .set_species("Rb")
                     .set_quantum_number_n(60)
                     .set_quantum_number_l(0)
                     .set_quantum_number_j(0.5)
                     .set_quantum_number_m(0.5)
                     .create(database);

    auto ket_p = KetAtomCreator<double>()
                     .set_species("Rb")
                     .set_quantum_number_n(60)
                     .set_quantum_number_l(1)
                     .set_quantum_number_j(0.5)
                     .set_quantum_number_m(0.5)
                     .create(database);

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(59, 61)
                     .restrict_quantum_number_l(0, 1)
                     .restrict_quantum_number_m(0.5, 0.5)
                     .create(database);

    SystemAtom<double> system(basis);

    DOCTEST_SUBCASE("calculate energy") {
        auto m1 = basis->get_canonical_state_from_ket(ket_s)->get_matrix_elements(
            ket_s, OperatorType::ENERGY, 0);
        DOCTEST_CHECK(m1.size() == 1);
        double energy1 = m1[0];

        auto m2 = basis->get_matrix_elements(ket_s, OperatorType::ENERGY, 0);
        DOCTEST_CHECK(m2.size() == basis->get_number_of_states());
        size_t idx = basis->get_corresponding_state_index(ket_s);
        double energy2 = m2[idx];

        double reference = ket_s->get_energy();
        DOCTEST_CHECK(std::abs(energy1 - reference) < 1e-11);
        DOCTEST_CHECK(std::abs(energy2 - reference) < 1e-11);
    }

    DOCTEST_SUBCASE("calculate electric dipole matrix element") {
        auto m = basis->get_matrix_elements(ket_p, OperatorType::ELECTRIC_DIPOLE, 0);
        DOCTEST_CHECK(m.size() == basis->get_number_of_states());
        size_t idx = basis->get_corresponding_state_index(ket_s);
        double dipole = m[idx];

        DOCTEST_CHECK(std::abs(dipole - 1247.5955234484584) < 1e-6);
    }

    DOCTEST_SUBCASE("calculate electric dipole matrix element with an induced dipole") {
        {
            auto state = basis->get_corresponding_state(ket_s);

            auto m = state->get_matrix_elements(state, OperatorType::ELECTRIC_DIPOLE, 0);
            DOCTEST_CHECK(m.rows() == 1);
            DOCTEST_CHECK(m.cols() == 1);
            double dipole = m.coeff(0, 0);

            DOCTEST_CHECK(std::abs(dipole - 0) < 1e-6);
        }

        {
            system.set_electric_field({0, 0, 1.9446903811524456e-10});
            system.diagonalize(DiagonalizerEigen<double>());
            auto state = system.get_eigenbasis()->get_corresponding_state(ket_s);

            auto m = state->get_matrix_elements(state, OperatorType::ELECTRIC_DIPOLE, 0);
            DOCTEST_CHECK(m.rows() == 1);
            DOCTEST_CHECK(m.cols() == 1);
            double dipole = m.coeff(0, 0);

            DOCTEST_CHECK(std::abs(dipole - 135.0169870551207) < 1e-6);
        }
    }
}

} // namespace pairinteraction
