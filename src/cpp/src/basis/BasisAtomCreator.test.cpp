// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtomCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <cmath>
#include <doctest/doctest.h>
#include <stdexcept>
#include <vector>

namespace pairinteraction {

constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;

DOCTEST_TEST_CASE("create a basis for strontium 88") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>()
                     .set_species("Sr88_sqdt")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 2)
                     .restrict_quantum_number("s", 0, 0)
                     .create(database);
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_sqdt");
    }
}

DOCTEST_TEST_CASE("create a basis for strontium 87") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>()
                     .set_species("Sr87_mqdt")
                     .restrict_quantum_number("nu", 59, 61)
                     .restrict_quantum_number("l", 0, 0)
                     .create(database);
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr87_mqdt");
    }
}

DOCTEST_TEST_CASE("create a basis from kets") {
    Database &database = Database::get_global_instance();
    auto ket1 = KetAtomCreator("Sr88_sqdt", 59, 0, 0, 0).create(database);
    auto ket2 = KetAtomCreator("Sr88_sqdt", 60, 0, 0, 0).create(database);
    auto ket3 = KetAtomCreator("Sr88_sqdt", 61, 0, 0, 0).create(database);
    auto basis =
        BasisAtomCreator<double>().add_ket(ket1).add_ket(ket2).add_ket(ket3).create(database);
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_sqdt");
    }
}

DOCTEST_TEST_CASE("create a basis and sort it according to parity and m") {
    Database &database = Database::get_global_instance();
    auto basis_unsorted = BasisAtomCreator<double>()
                              .set_species("Rb")
                              .restrict_quantum_number("n", 60, 60)
                              .restrict_quantum_number("l", 0, 3)
                              .restrict_quantum_number("m", -0.5, 0.5)
                              .create(database);

    // Sort the basis by parity and the m quantum number
    auto sorter = basis_unsorted->get_sorter({SorterType::PARITY, SorterType::QUANTUM_NUMBER_M});
    auto basis = basis_unsorted->transformed(sorter);

    // Check if the basis is properly sorted
    auto parity = Parity::ODD;
    auto quantum_number_m = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < basis->get_number_of_states(); ++i) {
        DOCTEST_MESSAGE("State ", i, ": Parity = ", basis->get_parity(i),
                        ", M = ", basis->get_quantum_number_m(i));
        DOCTEST_CHECK(basis->get_parity(i) >= parity);
        if (basis->get_parity(i) != parity) {
            parity = basis->get_parity(i);
            quantum_number_m = std::numeric_limits<double>::lowest();
        }
        DOCTEST_CHECK(basis->get_quantum_number_m(i) >= quantum_number_m);
        quantum_number_m = basis->get_quantum_number_m(i);
    }

    // Check that the blocks are correctly determined
    auto blocks = basis->get_indices_of_blocks({SorterType::PARITY, SorterType::QUANTUM_NUMBER_M});
    std::vector<size_t> expected_start = {0, 4, 8, 11};

    DOCTEST_CHECK(blocks.size() == expected_start.size());

    size_t idx = 0;
    for (const auto &block : blocks) {
        DOCTEST_MESSAGE("Block ", idx, " starts at ", block.start);
        DOCTEST_CHECK(block.start == expected_start[idx]);
        idx++;
    }

    // Transforming by the identity matrix leaves the coefficients unchanged
    size_t dim = basis->get_number_of_states();
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix(static_cast<long>(dim),
                                                        static_cast<long>(dim));
    matrix.setIdentity();
    auto transformed = basis->transformed(matrix);
    DOCTEST_CHECK(transformed->get_coefficients().isApprox(basis->get_coefficients()));
}

DOCTEST_TEST_CASE("a basis is canonical if its coefficients are the identity matrix") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 1)
                     .create(database);

    DOCTEST_CHECK(basis->is_canonical());

    auto dim = static_cast<long>(basis->get_number_of_states());

    // Explicitly stored zeros must not affect the result
    Eigen::SparseMatrix<double, Eigen::RowMajor> identity_with_zeros(dim, dim);
    identity_with_zeros.reserve(Eigen::VectorXi::Constant(dim, 2));
    for (long i = 0; i < dim; ++i) {
        identity_with_zeros.insert(i, i) = 1;
        identity_with_zeros.insert(i, (i + 1) % dim) = 0;
    }
    identity_with_zeros.makeCompressed();
    DOCTEST_CHECK(identity_with_zeros.nonZeros() > dim);
    DOCTEST_CHECK(basis->copy_with_coefficients(identity_with_zeros)->is_canonical());

    // A basis whose states are a non-trivial superposition of kets is not canonical
    Eigen::SparseMatrix<double, Eigen::RowMajor> non_identity(dim, dim);
    non_identity.setIdentity();
    non_identity.coeffRef(0, 0) = 0.5;
    DOCTEST_CHECK_FALSE(basis->copy_with_coefficients(non_identity)->is_canonical());

    // Sorting the basis makes it non-canonical
    auto sorted = basis->transformed(basis->get_sorter({SorterType::PARITY}));
    DOCTEST_CHECK_FALSE(sorted->is_canonical());
}

DOCTEST_TEST_CASE("blocks can only be obtained if the states are sorted") {
    Database &database = Database::get_global_instance();
    auto basis_unsorted = BasisAtomCreator<double>()
                              .set_species("Rb")
                              .restrict_quantum_number("n", 60, 60)
                              .restrict_quantum_number("l", 0, 3)
                              .restrict_quantum_number("m", -0.5, 0.5)
                              .create(database);

    // In the unsorted basis, states that share the same labels are scattered over the basis. They
    // would end up in several blocks and couplings between them would be lost.
    DOCTEST_CHECK_THROWS_AS(
        basis_unsorted->get_indices_of_blocks({SorterType::PARITY, SorterType::QUANTUM_NUMBER_M}),
        std::invalid_argument);

    auto basis = basis_unsorted->transformed(
        basis_unsorted->get_sorter({SorterType::PARITY, SorterType::QUANTUM_NUMBER_M}));

    // After sorting, the states of equal parity and m are contiguous ...
    DOCTEST_CHECK_NOTHROW(
        basis->get_indices_of_blocks({SorterType::PARITY, SorterType::QUANTUM_NUMBER_M}));

    // ... and so are the states of equal parity because the parity is the primary criterion
    DOCTEST_CHECK_NOTHROW(basis->get_indices_of_blocks({SorterType::PARITY}));

    // The states of equal m are not contiguous, however, because m only breaks ties between states
    // of equal parity
    DOCTEST_CHECK_THROWS_AS(basis->get_indices_of_blocks({SorterType::QUANTUM_NUMBER_M}),
                            std::invalid_argument);

    // If the states are not labeled at all, no blocks can be obtained either
    auto dim = static_cast<long>(basis->get_number_of_states());
    Eigen::SparseMatrix<double, Eigen::RowMajor> identity(dim, dim);
    identity.setIdentity();
    auto unlabeled = basis->copy_with_coefficients(identity);
    DOCTEST_CHECK_THROWS_AS(unlabeled->get_indices_of_blocks({SorterType::PARITY}),
                            std::invalid_argument);
    DOCTEST_CHECK_THROWS_AS(unlabeled->get_indices_of_blocks({SorterType::QUANTUM_NUMBER_F}),
                            std::invalid_argument);
}

DOCTEST_TEST_CASE("calculation of matrix elements") {
    auto &database = Database::get_global_instance();

    auto ket_s = KetAtomCreator()
                     .set_species("Rb")
                     .set_quantum_number("n", 60)
                     .set_quantum_number("l", 0)
                     .set_quantum_number("j", 0.5)
                     .set_quantum_number("m", 0.5)
                     .create(database);

    auto ket_p = KetAtomCreator()
                     .set_species("Rb")
                     .set_quantum_number("n", 60)
                     .set_quantum_number("l", 1)
                     .set_quantum_number("j", 0.5)
                     .set_quantum_number("m", 0.5)
                     .create(database);

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number("n", 59, 61)
                     .restrict_quantum_number("l", 0, 1)
                     .restrict_quantum_number("m", 0.5, 0.5)
                     .create(database);

    SystemAtom<double> system(basis);

    auto get_corresponding_state_index = [&database](const auto &b,
                                                     const std::shared_ptr<const KetAtom> &ket) {
        auto basis_ket = BasisAtomCreator<double>().add_ket(ket).create(database);
        Eigen::MatrixXd overlaps =
            Eigen::MatrixXd(b->get_matrix_elements(basis_ket, OperatorType::IDENTITY, 0))
                .cwiseAbs();
        Eigen::Index idx = 0;
        overlaps.row(0).maxCoeff(&idx);
        return idx;
    };

    DOCTEST_SUBCASE("calculate energy") {
        auto basis_ket_s = BasisAtomCreator<double>().add_ket(ket_s).create(database);

        auto m1 = basis_ket_s->get_matrix_elements(basis_ket_s, OperatorType::ENERGY, 0);
        DOCTEST_CHECK(m1.rows() == 1);
        DOCTEST_CHECK(m1.cols() == 1);
        double energy1 = m1.coeff(0, 0);

        auto m2 = basis->get_matrix_elements(basis_ket_s, OperatorType::ENERGY, 0);
        DOCTEST_CHECK(m2.rows() == 1);
        DOCTEST_CHECK(m2.cols() == basis->get_number_of_states());
        double energy2 = m2.coeff(0, static_cast<int>(get_corresponding_state_index(basis, ket_s)));

        double reference = ket_s->get_energy();
        DOCTEST_CHECK(std::abs(energy1 - reference) < 1e-11);
        DOCTEST_CHECK(std::abs(energy2 - reference) < 1e-11);
    }

    DOCTEST_SUBCASE("calculate electric dipole matrix element") {
        auto basis_ket_p = BasisAtomCreator<double>().add_ket(ket_p).create(database);

        auto m = basis->get_matrix_elements(basis_ket_p, OperatorType::ELECTRIC_DIPOLE, 0);
        DOCTEST_CHECK(m.rows() == 1);
        DOCTEST_CHECK(m.cols() == basis->get_number_of_states());
        double dipole = m.coeff(0, static_cast<int>(get_corresponding_state_index(basis, ket_s)));

        DOCTEST_CHECK(std::abs(dipole - 1247.60438) < 1e-3);
    }

    DOCTEST_SUBCASE("calculate electric dipole matrix element with and without an induced dipole") {
        {
            auto state = basis->get_state(get_corresponding_state_index(basis, ket_s));

            auto m = state->get_matrix_elements(state, OperatorType::ELECTRIC_DIPOLE, 0);
            DOCTEST_CHECK(m.rows() == 1);
            DOCTEST_CHECK(m.cols() == 1);
            double dipole = m.coeff(0, 0);

            DOCTEST_CHECK(std::abs(dipole - 0) < 1e-6);
        }

        {
            system.set_electric_field({0, 0, VOLT_PER_CM_IN_ATOMIC_UNITS});
            system.diagonalize(DiagonalizerEigen<double>());
            auto eigenbasis = system.get_eigenbasis();
            auto state = eigenbasis->get_state(get_corresponding_state_index(eigenbasis, ket_s));

            auto m = state->get_matrix_elements(state, OperatorType::ELECTRIC_DIPOLE, 0);
            DOCTEST_CHECK(m.rows() == 1);
            DOCTEST_CHECK(m.cols() == 1);
            double dipole = m.coeff(0, 0);

            DOCTEST_CHECK(std::abs(dipole - 135.04130) < 1e-3);
        }
    }
}

DOCTEST_TEST_CASE("conserved quantum numbers are detected despite non-normalized coefficients") {
    // The diagonalizers prune small entries of the eigenvectors so that the columns of the
    // resulting transformation matrix are normalized only up to a defect of the order of rtol^2.
    // The detection of conserved quantum numbers must not be spoiled by this defect. Otherwise,
    // the detection would fail for a large rtol or a low floating point precision.
    Database &database = Database::get_global_instance();
    auto unsorted_basis = BasisAtomCreator<double>()
                              .set_species("Rb")
                              .restrict_quantum_number("n", 60, 60)
                              .restrict_quantum_number("l", 0, 3)
                              .create(database);
    auto basis =
        unsorted_basis->transformed(unsorted_basis->get_sorter({SorterType::QUANTUM_NUMBER_M}));
    DOCTEST_REQUIRE(basis->has_quantum_number_m());

    auto dim = static_cast<Eigen::Index>(basis->get_number_of_states());
    constexpr double defect = 1e-2; // the squared norm of each column is 1 - defect
    const double scale = std::sqrt(1 - defect);
    const double sqrt_half = std::sqrt(0.5);

    // Superimpose pairs of states that have the same m, mimicking eigenvectors that are confined
    // to a block of constant m, and scale the columns so that they are not exactly normalized
    std::vector<Eigen::Triplet<double>> triplets;
    for (Eigen::Index i = 0; i < dim;) {
        auto idx = static_cast<size_t>(i);
        if (i + 1 < dim &&
            basis->get_quantum_number_m(idx) == basis->get_quantum_number_m(idx + 1)) {
            triplets.emplace_back(i, i, scale * sqrt_half);
            triplets.emplace_back(i + 1, i, scale * sqrt_half);
            triplets.emplace_back(i, i + 1, scale * sqrt_half);
            triplets.emplace_back(i + 1, i + 1, -scale * sqrt_half);
            i += 2;
        } else {
            triplets.emplace_back(i, i, scale);
            i += 1;
        }
    }
    Eigen::SparseMatrix<double, Eigen::RowMajor> transformation(dim, dim);
    transformation.setFromTriplets(triplets.begin(), triplets.end());

    // Because the transformation does not mix different m, m is still well-defined
    auto transformed = basis->transformed(transformation);
    DOCTEST_CHECK(transformed->has_quantum_number_m());
    for (size_t i = 0; i < transformed->get_number_of_states(); ++i) {
        DOCTEST_CHECK(transformed->get_quantum_number_m(i) == basis->get_quantum_number_m(i));
    }

    // If, in contrast, states of different m are superimposed, m is not well-defined anymore
    Eigen::Index other = 0;
    while (other < dim &&
           basis->get_quantum_number_m(0) ==
               basis->get_quantum_number_m(static_cast<size_t>(other))) {
        ++other;
    }
    DOCTEST_REQUIRE(other < dim);
    triplets.clear();
    for (Eigen::Index i = 0; i < dim; ++i) {
        if (i != 0 && i != other) {
            triplets.emplace_back(i, i, scale);
        }
    }
    triplets.emplace_back(0, 0, scale * sqrt_half);
    triplets.emplace_back(other, 0, scale * sqrt_half);
    triplets.emplace_back(0, other, scale * sqrt_half);
    triplets.emplace_back(other, other, -scale * sqrt_half);
    Eigen::SparseMatrix<double, Eigen::RowMajor> mixing_transformation(dim, dim);
    mixing_transformation.setFromTriplets(triplets.begin(), triplets.end());

    auto mixed = basis->transformed(mixing_transformation);
    DOCTEST_CHECK_FALSE(mixed->has_quantum_number_m());
    DOCTEST_CHECK_THROWS_AS(mixed->get_quantum_number_m(0), std::invalid_argument);
}

} // namespace pairinteraction
