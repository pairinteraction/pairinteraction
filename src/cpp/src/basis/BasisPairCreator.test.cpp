// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPairCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetAtomCreator.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/system/SystemPair.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <doctest/doctest.h>

namespace pairinteraction {

constexpr double HARTREE_IN_GHZ = 6579683.920501762;
constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;

DOCTEST_TEST_CASE("create a BasisPair") {
    // Create single-atom system
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    SystemAtom<double> system(basis);
    system.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});

    DiagonalizerEigen<double> diagonalizer;
    system.diagonalize(diagonalizer);

    // Get energy window for a two-atom basis
    auto ket = KetAtomCreator()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);
    double min_energy = 2 * ket->get_energy() - 3 / HARTREE_IN_GHZ;
    double max_energy = 2 * ket->get_energy() + 3 / HARTREE_IN_GHZ;

    // Create two-atom bases
    auto basis_pair_a = pairinteraction::BasisPairCreator<double>()
                            .add(system)
                            .add(system)
                            .restrict_energy(min_energy, max_energy)
                            .restrict_quantum_number_m(1, 1)
                            .create();
    auto basis_pair_b = pairinteraction::BasisPairCreator<double>()
                            .add(system)
                            .add(system)
                            .restrict_energy(min_energy, max_energy)
                            .restrict_quantum_number_m(1, 1)
                            .create();

    DOCTEST_SUBCASE("check equality of kets") {
        // Obtain kets from the two-atom bases and check for equality
        auto ket1a = basis_pair_a->get_kets()[0];
        auto ket1b = basis_pair_b->get_kets()[0];
        auto ket2a = basis_pair_a->get_kets()[1];
        auto ket2b = basis_pair_b->get_kets()[1];
        DOCTEST_CHECK(*ket1a == *ket1a);
        DOCTEST_CHECK(*ket2a == *ket2a);
        DOCTEST_CHECK(*ket1a != *ket2b);
        DOCTEST_CHECK(*ket2a != *ket1b);

        // Currently, kets from different BasisPair are never equal
        DOCTEST_CHECK(*ket1a != *ket1b);
        DOCTEST_CHECK(*ket2a != *ket2b);
    }

    DOCTEST_SUBCASE("check overlap") {
        auto overlaps = basis_pair_a->get_overlaps(ket, ket);

        // The total overlap is less than 1 because of the restricted energy window
        DOCTEST_CHECK(overlaps.sum() == doctest::Approx(0.9107819201));
    }

    DOCTEST_SUBCASE("get the atomic states constituting a ket of the basis_pair") {
        auto atomic_states = basis_pair_a->get_kets()[0]->get_atomic_states();
        DOCTEST_CHECK(atomic_states.size() == 2);
        DOCTEST_CHECK(atomic_states[0]->get_number_of_states() == 1);
        DOCTEST_CHECK(atomic_states[0]->get_number_of_kets() == basis->get_number_of_kets());
    }
}

DOCTEST_TEST_CASE("get matrix elements in the pair basis") {
    DiagonalizerEigen<double> diagonalizer;

    // Create single-atom system
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 2)
                     .create(database);
    SystemAtom<double> system(basis);
    system.set_electric_field({0, 0, 10 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    system.diagonalize(diagonalizer);

    // Get energy window for a two-atom basis
    auto ket = KetAtomCreator()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);
    double min_energy = 2 * ket->get_energy() - 3 / HARTREE_IN_GHZ;
    double max_energy = 2 * ket->get_energy() + 3 / HARTREE_IN_GHZ;

    // Create two-atom system
    auto basis_pair_unperturbed = pairinteraction::BasisPairCreator<double>()
                                      .add(system)
                                      .add(system)
                                      .restrict_energy(min_energy, max_energy)
                                      .restrict_quantum_number_m(1, 1)
                                      .create();
    auto system_pair = SystemPair<double>(basis_pair_unperturbed)
                           .set_distance_vector({0, 0, 1 * UM_IN_ATOMIC_UNITS});
    system_pair.diagonalize(diagonalizer);

    auto basis_pair = system_pair.get_eigenbasis();

    DOCTEST_SUBCASE("check dimensions") {
        // <basis_pair_unperturbed|d0d0|basis_pair_unperturbed>
        auto matrix_elements_all = basis_pair_unperturbed->get_matrix_elements(
            basis_pair_unperturbed, OperatorType::ELECTRIC_DIPOLE, OperatorType::ELECTRIC_DIPOLE, 0,
            0);
        DOCTEST_CHECK(matrix_elements_all.rows() == basis_pair_unperturbed->get_number_of_states());
        DOCTEST_CHECK(matrix_elements_all.cols() == basis_pair_unperturbed->get_number_of_states());

        // <ket_pair|d0d0|basis_pair_unperturbed>
        auto ket_pair = basis_pair_unperturbed->get_kets()[0];
        auto matrix_elements_ket_pair = basis_pair_unperturbed->get_matrix_elements(
            ket_pair, OperatorType::ELECTRIC_DIPOLE, OperatorType::ELECTRIC_DIPOLE, 0, 0);
        DOCTEST_CHECK(matrix_elements_ket_pair.size() ==
                      basis_pair_unperturbed->get_number_of_states());

        {
            Eigen::VectorX<double> ref = matrix_elements_all.row(0);
            DOCTEST_CHECK(ref.isApprox(matrix_elements_ket_pair, 1e-11));
        }

        // <basis x basis|d0d0|basis_pair>
        auto matrix_elements_product = basis_pair->get_matrix_elements(
            basis, basis, OperatorType::ELECTRIC_DIPOLE, OperatorType::ELECTRIC_DIPOLE, 0, 0);
        DOCTEST_CHECK(matrix_elements_product.rows() ==
                      basis->get_number_of_states() * basis->get_number_of_states());
        DOCTEST_CHECK(matrix_elements_product.cols() == basis_pair->get_number_of_states());

        // <ket,ket|d0d0|basis_pair>
        auto matrix_elements_ket = basis_pair->get_matrix_elements(
            ket, ket, OperatorType::ELECTRIC_DIPOLE, OperatorType::ELECTRIC_DIPOLE, 0, 0);
        DOCTEST_CHECK(matrix_elements_ket.size() == basis_pair->get_number_of_states());
    }

    DOCTEST_SUBCASE("check matrix elements") {
        // energy
        auto hamiltonian = basis_pair->get_matrix_elements(basis_pair, OperatorType::ENERGY,
                                                           OperatorType::IDENTITY);
        hamiltonian += basis_pair->get_matrix_elements(basis_pair, OperatorType::IDENTITY,
                                                       OperatorType::ENERGY);

        // interaction with electric field
        {
            Eigen::SparseMatrix<double, Eigen::RowMajor> tmp = -basis_pair->get_matrix_elements(
                basis_pair, OperatorType::ELECTRIC_DIPOLE, OperatorType::IDENTITY, 0, 0);
            tmp += -basis_pair->get_matrix_elements(basis_pair, OperatorType::IDENTITY,
                                                    OperatorType::ELECTRIC_DIPOLE, 0, 0);
            hamiltonian += 10 * VOLT_PER_CM_IN_ATOMIC_UNITS * tmp;
        }

        // dipole-dipole interaction
        {
            Eigen::SparseMatrix<double, Eigen::RowMajor> tmp = -2 *
                basis_pair->get_matrix_elements(basis_pair, OperatorType::ELECTRIC_DIPOLE,
                                                OperatorType::ELECTRIC_DIPOLE, 0, 0);
            tmp += -basis_pair->get_matrix_elements(basis_pair, OperatorType::ELECTRIC_DIPOLE,
                                                    OperatorType::ELECTRIC_DIPOLE, 1, -1);
            tmp += -basis_pair->get_matrix_elements(basis_pair, OperatorType::ELECTRIC_DIPOLE,
                                                    OperatorType::ELECTRIC_DIPOLE, -1, 1);
            hamiltonian += std::pow(UM_IN_ATOMIC_UNITS, -3) * tmp;
        }

        // compare to reference
        const auto &ref = system_pair.get_matrix();
        DOCTEST_CHECK(ref.isApprox(hamiltonian, 1e-11));
    }
}

} // namespace pairinteraction
