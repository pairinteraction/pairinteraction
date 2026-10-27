// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisAtomCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/StateLabel.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <cmath>
#include <doctest/doctest.h>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
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
    auto sorter = basis_unsorted->get_sorter({StateLabel::PARITY, StateLabel::QUANTUM_NUMBER_M});
    auto basis = basis_unsorted->get_sorted_copy(sorter);

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
    auto blocks = basis->get_indices_of_blocks({StateLabel::PARITY, StateLabel::QUANTUM_NUMBER_M});
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
    DOCTEST_CHECK(basis->get_modified_copy(identity_with_zeros, std::nullopt, {})->is_canonical());

    // A basis whose states are a non-trivial superposition of kets is not canonical
    Eigen::SparseMatrix<double, Eigen::RowMajor> non_identity(dim, dim);
    non_identity.setIdentity();
    non_identity.coeffRef(0, 0) = 0.5;
    DOCTEST_CHECK_FALSE(basis->get_modified_copy(non_identity, std::nullopt, {})->is_canonical());

    // Sorting the basis makes it non-canonical
    auto sorted = basis->get_sorted_copy(basis->get_sorter({StateLabel::PARITY}));
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
        basis_unsorted->get_indices_of_blocks({StateLabel::PARITY, StateLabel::QUANTUM_NUMBER_M}),
        std::invalid_argument);

    auto basis = basis_unsorted->get_sorted_copy(
        basis_unsorted->get_sorter({StateLabel::PARITY, StateLabel::QUANTUM_NUMBER_M}));

    // After sorting, the states of equal parity and m are contiguous ...
    DOCTEST_CHECK_NOTHROW(
        basis->get_indices_of_blocks({StateLabel::PARITY, StateLabel::QUANTUM_NUMBER_M}));

    // ... and so are the states of equal parity because the parity is the primary criterion
    DOCTEST_CHECK_NOTHROW(basis->get_indices_of_blocks({StateLabel::PARITY}));

    // The states of equal m are not contiguous, however, because m only breaks ties between states
    // of equal parity
    DOCTEST_CHECK_THROWS_AS(basis->get_indices_of_blocks({StateLabel::QUANTUM_NUMBER_M}),
                            std::invalid_argument);

    // If the states are not labeled at all, no blocks can be obtained either
    auto dim = static_cast<long>(basis->get_number_of_states());
    Eigen::SparseMatrix<double, Eigen::RowMajor> identity(dim, dim);
    identity.setIdentity();
    auto unlabeled = basis->get_modified_copy(identity, std::nullopt, {});
    DOCTEST_CHECK_THROWS_AS(unlabeled->get_indices_of_blocks({StateLabel::PARITY}),
                            std::invalid_argument);
    DOCTEST_CHECK_THROWS_AS(unlabeled->get_indices_of_blocks({StateLabel::QUANTUM_NUMBER_F}),
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
    auto sorter = basis->get_sorter({StateLabel::QUANTUM_NUMBER_M});
    auto sorted = basis->get_sorted_copy(sorter);
    DOCTEST_CHECK(sorted->has_energy());
    DOCTEST_CHECK_FALSE(sorted->is_sorted_by_energy());
    for (Eigen::Index i = 0; i < dim; ++i) {
        DOCTEST_CHECK(sorted->get_energy(static_cast<size_t>(i)) ==
                      basis->get_energy(static_cast<size_t>(sorter.indices()[i])));
    }

    // Sorting by the energy restores the order of the energies
    auto resorted = sorted->get_sorted_copy(sorted->get_sorter_by_energy());
    DOCTEST_CHECK(resorted->is_sorted_by_energy());

    // A modified copy of a basis has the explicitly specified energies
    std::vector<double> energies(basis->get_number_of_states());
    for (size_t i = 0; i < energies.size(); ++i) {
        energies[i] = basis->get_energy(i);
    }
    auto with_energies = basis->get_modified_copy(basis->get_coefficients(), energies,
                                                  basis->get_quantum_numbers_of_states());
    DOCTEST_CHECK(with_energies->is_sorted_by_energy());
    for (size_t i = 0; i < energies.size(); ++i) {
        DOCTEST_CHECK(with_energies->get_energy(i) == energies[i]);
    }
    DOCTEST_CHECK_THROWS_AS(basis->get_modified_copy(basis->get_coefficients(),
                                                     std::vector<double>{energies[0]},
                                                     basis->get_quantum_numbers_of_states()),
                            std::invalid_argument);

    // Energies that are not given are not well-defined
    auto copied = basis->get_modified_copy(basis->get_coefficients(), std::nullopt,
                                           basis->get_quantum_numbers_of_states());
    DOCTEST_CHECK_THROWS_AS(copied->get_energy(0), std::invalid_argument);
    DOCTEST_CHECK_FALSE(copied->has_energy());
    DOCTEST_CHECK_FALSE(copied->is_sorted_by_energy());
    DOCTEST_CHECK_THROWS_AS(copied->get_sorter_by_energy(), std::invalid_argument);
    DOCTEST_CHECK_FALSE(copied->get_state(0)->has_energy());

    // Energies that are given must be finite and smaller than the maximum value
    for (double invalid_energy :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()}) {
        energies[0] = invalid_energy;
        DOCTEST_CHECK_THROWS_AS(basis->get_modified_copy(basis->get_coefficients(), energies,
                                                         basis->get_quantum_numbers_of_states()),
                                std::invalid_argument);
    }
}

DOCTEST_TEST_CASE("quantum numbers of the states of a basis are (half-)integers") {
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 60)
                     .restrict_quantum_number("l", 0, 1)
                     .create();

    auto with_m = [&](double first_m) {
        auto quantum_numbers_of_states = basis->get_quantum_numbers_of_states();
        quantum_numbers_of_states.at("m")[0] = first_m;
        return basis->get_modified_copy(basis->get_coefficients(), std::nullopt,
                                        std::move(quantum_numbers_of_states));
    };

    // Integers and half-integers are allowed
    DOCTEST_CHECK_NOTHROW(with_m(-1.5));
    DOCTEST_CHECK_NOTHROW(with_m(2));

    // Other values cannot be doubled to an exact integer. In particular, there is no value that
    // marks a quantum number as not well-defined for a single state.
    DOCTEST_CHECK_THROWS_AS(with_m(0.3), std::invalid_argument);
    DOCTEST_CHECK_THROWS_AS(with_m(0.5 + 1e-12), std::invalid_argument);
    DOCTEST_CHECK_THROWS_AS(with_m(std::numeric_limits<double>::quiet_NaN()),
                            std::invalid_argument);
    DOCTEST_CHECK_THROWS_AS(with_m(std::numeric_limits<double>::infinity()), std::invalid_argument);
    DOCTEST_CHECK_THROWS_AS(with_m(std::numeric_limits<double>::max()), std::invalid_argument);
    DOCTEST_CHECK_THROWS_AS(with_m(std::numeric_limits<double>::lowest()), std::invalid_argument);

    // A quantum number must be given for all states
    auto quantum_numbers_of_states = basis->get_quantum_numbers_of_states();
    quantum_numbers_of_states.at("m").pop_back();
    DOCTEST_CHECK_THROWS_AS(basis->get_modified_copy(basis->get_coefficients(), std::nullopt,
                                                     std::move(quantum_numbers_of_states)),
                            std::invalid_argument);
}

} // namespace pairinteraction
