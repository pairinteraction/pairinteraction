// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtomCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
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
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Sr88_sqdt")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 2)
                     .restrict_quantum_number("s", 0, 0)
                     .create();
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_sqdt");
    }
}

DOCTEST_TEST_CASE("create a basis for strontium 87") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Sr87_mqdt")
                     .restrict_quantum_number("nu", 59, 61)
                     .restrict_quantum_number("l", 0, 0)
                     .create();
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr87_mqdt");
    }
}

DOCTEST_TEST_CASE("create a basis from kets") {
    Database &database = Database::get_global_instance();
    auto ket1 = KetAtomCreator(database)
                    .set_species("Sr88_sqdt")
                    .set_quantum_number("n", 59)
                    .set_quantum_number("l", 0)
                    .set_quantum_number("j", 0)
                    .set_quantum_number("m", 0)
                    .create();
    auto ket2 = KetAtomCreator(database)
                    .set_species("Sr88_sqdt")
                    .set_quantum_number("n", 60)
                    .set_quantum_number("l", 0)
                    .set_quantum_number("j", 0)
                    .set_quantum_number("m", 0)
                    .create();
    auto ket3 = KetAtomCreator(database)
                    .set_species("Sr88_sqdt")
                    .set_quantum_number("n", 61)
                    .set_quantum_number("l", 0)
                    .set_quantum_number("j", 0)
                    .set_quantum_number("m", 0)
                    .create();
    auto basis =
        BasisAtomCreator<double>(database).add_ket(ket1).add_ket(ket2).add_ket(ket3).create();
    for (const auto &ket : *basis) {
        DOCTEST_CHECK(ket->get_species() == "Sr88_sqdt");
    }
}

DOCTEST_TEST_CASE("create a basis and sort it according to parity and m") {
    Database &database = Database::get_global_instance();
    auto basis_unsorted = BasisAtomCreator<double>(database)
                              .set_species("Rb")
                              .restrict_quantum_number("n", 60, 60)
                              .restrict_quantum_number("l", 0, 3)
                              .restrict_quantum_number("m", -0.5, 0.5)
                              .create();

    // Sort the basis by parity and the m quantum number
    auto sorter = basis_unsorted->get_sorter({SorterType::PARITY, SorterType::QUANTUM_NUMBER_M});
    auto basis = basis_unsorted->sorted(sorter);

    // Check if the basis is properly sorted
    double parity = -1;
    auto quantum_number_m = std::numeric_limits<double>::lowest();
    for (size_t i = 0; i < basis->get_number_of_states(); ++i) {
        DOCTEST_MESSAGE("State ", i, ": Parity = ", basis->get_quantum_number("parity", i),
                        ", M = ", basis->get_quantum_number("m", i));
        DOCTEST_CHECK(basis->get_quantum_number("parity", i) >= parity);
        if (basis->get_quantum_number("parity", i) != parity) {
            parity = basis->get_quantum_number("parity", i);
            quantum_number_m = std::numeric_limits<double>::lowest();
        }
        DOCTEST_CHECK(basis->get_quantum_number("m", i) >= quantum_number_m);
        quantum_number_m = basis->get_quantum_number("m", i);
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
}

DOCTEST_TEST_CASE("a basis is canonical if its coefficients are the identity matrix") {
    auto &database = Database::get_global_instance();

    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 1)
                     .create();

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
    auto sorted = basis->sorted(basis->get_sorter({SorterType::PARITY}));
    DOCTEST_CHECK_FALSE(sorted->is_canonical());
}

DOCTEST_TEST_CASE("blocks can only be obtained if the states are sorted") {
    Database &database = Database::get_global_instance();
    auto basis_unsorted = BasisAtomCreator<double>(database)
                              .set_species("Rb")
                              .restrict_quantum_number("n", 60, 60)
                              .restrict_quantum_number("l", 0, 3)
                              .restrict_quantum_number("m", -0.5, 0.5)
                              .create();

    // In the unsorted basis, states that share the same labels are scattered over the basis. They
    // would end up in several blocks and couplings between them would be lost.
    DOCTEST_CHECK_THROWS_AS(
        basis_unsorted->get_indices_of_blocks({SorterType::PARITY, SorterType::QUANTUM_NUMBER_M}),
        std::invalid_argument);

    auto basis = basis_unsorted->sorted(
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

    auto ket_s = KetAtomCreator(database)
                     .set_species("Rb")
                     .set_quantum_number("n", 60)
                     .set_quantum_number("l", 0)
                     .set_quantum_number("j", 0.5)
                     .set_quantum_number("m", 0.5)
                     .create();

    auto ket_p = KetAtomCreator(database)
                     .set_species("Rb")
                     .set_quantum_number("n", 60)
                     .set_quantum_number("l", 1)
                     .set_quantum_number("j", 0.5)
                     .set_quantum_number("m", 0.5)
                     .create();

    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 59, 61)
                     .restrict_quantum_number("l", 0, 1)
                     .restrict_quantum_number("m", 0.5, 0.5)
                     .create();

    SystemAtom<double> system(basis);

    auto get_corresponding_state_index = [&database](const auto &b,
                                                     const std::shared_ptr<const KetAtom> &ket) {
        auto basis_ket = BasisAtomCreator<double>(database).add_ket(ket).create();
        Eigen::MatrixXd overlaps =
            Eigen::MatrixXd(b->get_matrix_elements(basis_ket, OperatorType::IDENTITY, 0))
                .cwiseAbs();
        Eigen::Index idx = 0;
        overlaps.row(0).maxCoeff(&idx);
        return idx;
    };

    DOCTEST_SUBCASE("calculate energy") {
        auto basis_ket_s = BasisAtomCreator<double>(database).add_ket(ket_s).create();

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
        auto basis_ket_p = BasisAtomCreator<double>(database).add_ket(ket_p).create();

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

DOCTEST_TEST_CASE("energies of the states of a basis") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 1)
                     .create();
    auto dim = static_cast<Eigen::Index>(basis->get_number_of_states());

    // The states of a canonical basis have the energies of the kets, by which they are sorted
    DOCTEST_CHECK(basis->has_energy());
    DOCTEST_CHECK(basis->is_sorted_by_energy());
    for (size_t i = 0; i < basis->get_number_of_states(); ++i) {
        DOCTEST_CHECK(basis->get_energy(i) == basis->get_ket(i)->get_energy());
        DOCTEST_CHECK(basis->get_state(i)->get_energy(0) == basis->get_ket(i)->get_energy());
    }

    // Sorting permutes the energies along with the states
    auto sorter = basis->get_sorter({SorterType::QUANTUM_NUMBER_M});
    auto sorted = basis->sorted(sorter);
    DOCTEST_CHECK(sorted->has_energy());
    DOCTEST_CHECK_FALSE(sorted->is_sorted_by_energy());
    for (Eigen::Index i = 0; i < dim; ++i) {
        DOCTEST_CHECK(sorted->get_energy(static_cast<size_t>(i)) ==
                      basis->get_energy(static_cast<size_t>(sorter.indices()[i])));
    }

    // States that are constructed with unknown energies do not have well-defined energies
    std::vector<double> unknown_energies(static_cast<size_t>(dim),
                                         std::numeric_limits<double>::max());
    auto basis_without_energies = basis->copy_with_states(
        basis->get_coefficients(), unknown_energies, basis->get_quantum_numbers());
    DOCTEST_CHECK_FALSE(basis_without_energies->has_energy());
    DOCTEST_CHECK_FALSE(basis_without_energies->is_sorted_by_energy());
    DOCTEST_CHECK_THROWS_AS(basis_without_energies->get_energy(0), std::invalid_argument);

    // Otherwise, the states have the energies with which they are constructed
    std::vector<double> energies(basis->get_number_of_states());
    for (size_t i = 0; i < energies.size(); ++i) {
        energies[i] = 2 * basis->get_energy(i);
    }
    auto basis_with_energies =
        basis->copy_with_states(basis->get_coefficients(), energies, basis->get_quantum_numbers());
    DOCTEST_CHECK(basis_with_energies->is_sorted_by_energy());
    for (size_t i = 0; i < energies.size(); ++i) {
        DOCTEST_CHECK(basis_with_energies->get_energy(i) == energies[i]);
    }

    // Arbitrary coefficients do not have a well-defined energy
    auto copied = basis->copy_with_coefficients(basis->get_coefficients());
    DOCTEST_CHECK_THROWS_AS(copied->get_energy(0), std::invalid_argument);
    DOCTEST_CHECK_FALSE(copied->has_energy());
}

DOCTEST_TEST_CASE("copy a basis with given coefficients, energies, and quantum numbers") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 1)
                     .create();
    auto num_kets = static_cast<Eigen::Index>(basis->get_number_of_kets());

    // Construct a basis whose only state is a superposition of the first two kets
    Eigen::SparseMatrix<double, Eigen::RowMajor> coefficients(num_kets, 1);
    coefficients.insert(0, 0) = 1 / std::sqrt(2.0);
    coefficients.insert(1, 0) = 1 / std::sqrt(2.0);
    std::vector<double> energies{0.5};
    // Quantum numbers that are not given, like f, are not well-defined
    typename BasisAtom<double>::quantum_numbers_t quantum_numbers{{"m", {0.5}}, {"parity", {-1}}};

    auto constructed = basis->copy_with_states(coefficients, energies, quantum_numbers);
    DOCTEST_CHECK(constructed->get_number_of_kets() == basis->get_number_of_kets());
    DOCTEST_CHECK(constructed->get_number_of_states() == 1);
    DOCTEST_CHECK(constructed->get_coefficients().isApprox(coefficients));
    DOCTEST_CHECK(constructed->get_energy(0) == 0.5);
    DOCTEST_CHECK(constructed->get_quantum_number("m", 0) == 0.5);
    DOCTEST_CHECK(constructed->get_quantum_number("parity", 0) == -1);
    DOCTEST_CHECK_FALSE(constructed->has_quantum_number("f"));
    DOCTEST_CHECK(constructed->get_species() == basis->get_species());

    // The number of rows of the coefficients must match the number of kets
    Eigen::SparseMatrix<double, Eigen::RowMajor> too_few_rows(num_kets - 1, 1);
    DOCTEST_CHECK_THROWS_AS(basis->copy_with_states(too_few_rows, energies, quantum_numbers),
                            std::invalid_argument);

    // There must be one energy per state
    DOCTEST_CHECK_THROWS_AS(basis->copy_with_states(coefficients, {0.5, 0.5}, quantum_numbers),
                            std::invalid_argument);

    // The states can only be labeled by quantum numbers by which they can be sorted, ...
    auto unknown_quantum_number = quantum_numbers;
    unknown_quantum_number["n"] = {60};
    DOCTEST_CHECK_THROWS_AS(basis->copy_with_states(coefficients, energies, unknown_quantum_number),
                            std::invalid_argument);
    // ... and each given quantum number must have one value per state
    auto wrong_size = quantum_numbers;
    wrong_size["m"] = {0.5, 0.5};
    DOCTEST_CHECK_THROWS_AS(basis->copy_with_states(coefficients, energies, wrong_size),
                            std::invalid_argument);
}

} // namespace pairinteraction
