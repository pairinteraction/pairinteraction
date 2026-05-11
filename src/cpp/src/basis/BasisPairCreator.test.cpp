// SPDX-FileCopyrightText: 2024 PairInteraction Developers
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
#include "pairinteraction/utils/hash.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <doctest/doctest.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pairinteraction {

constexpr double HARTREE_IN_GHZ = 6579683.920501762;
constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;

namespace {
template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
build_manual_symmetrizer(const std::shared_ptr<const BasisPair<Scalar>> &basis,
                         Parity parity_under_inversion, Parity parity_under_permutation) {
    using real_t = typename BasisPair<Scalar>::real_t;

    const auto basis1 = basis->get_basis1();
    const auto basis2 = basis->get_basis2();
    const auto inv_sqrt_two = static_cast<real_t>(1 / std::sqrt(2.0));

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(2 * basis->get_number_of_states());

    Eigen::Index state_index = 0;
    std::unordered_map<std::array<size_t, 2>, Eigen::Index, utils::hash<std::array<size_t, 2>>>
        ket_ids_to_state_index;
    for (size_t idx1 = 0; idx1 < basis1->get_number_of_states(); ++idx1) {
        for (size_t idx2 = 0; idx2 < basis2->get_number_of_states(); ++idx2) {
            int ket_index = basis->get_ket_index_from_tuple(idx1, idx2);
            if (ket_index < 0) {
                continue;
            }

            size_t id1 = basis1->get_corresponding_ket(idx1)->get_id_in_database();
            size_t id2 = basis2->get_corresponding_ket(idx2)->get_id_in_database();

            if (id1 == id2) {
                if (parity_under_inversion == Parity::EVEN ||
                    parity_under_permutation == Parity::EVEN) {
                    continue;
                }
                triplets.emplace_back(ket_index, state_index++, Scalar{1});
                continue;
            }

            std::array<size_t, 2> ordered_ids{std::max(id1, id2), std::min(id1, id2)};
            auto [it, inserted] = ket_ids_to_state_index.try_emplace(ordered_ids, state_index);
            if (inserted) {
                ++state_index;
            }
            Eigen::Index column_index = it->second;

            if (id1 > id2) {
                triplets.emplace_back(ket_index, column_index, static_cast<Scalar>(inv_sqrt_two));
            } else {
                int swapped_sign = 0;
                if (parity_under_inversion != Parity::UNKNOWN &&
                    parity_under_permutation == Parity::UNKNOWN) {
                    swapped_sign = -static_cast<int>(parity_under_inversion) *
                        static_cast<int>(basis1->get_parity(idx1)) *
                        static_cast<int>(basis2->get_parity(idx2));
                } else {
                    swapped_sign = -static_cast<int>(parity_under_permutation);
                }
                triplets.emplace_back(ket_index, column_index,
                                      static_cast<Scalar>(swapped_sign * inv_sqrt_two));
            }
        }
    }

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> transformation(
        static_cast<Eigen::Index>(basis->get_number_of_states()), state_index);
    transformation.setFromTriplets(triplets.begin(), triplets.end());
    return transformation;
}

template <typename Scalar>
void check_same_pair_eigenenergies(const std::shared_ptr<const BasisPair<Scalar>> &basis1,
                                   const std::shared_ptr<const BasisPair<Scalar>> &basis2,
                                   const DiagonalizerEigen<Scalar> &diagonalizer) {
    auto system_pair_1 = SystemPair<Scalar>(basis1).set_distance_vector(
        std::array<typename BasisPair<Scalar>::real_t, 3>{0, 0, 1 * UM_IN_ATOMIC_UNITS});
    auto system_pair_2 = SystemPair<Scalar>(basis2).set_distance_vector(
        std::array<typename BasisPair<Scalar>::real_t, 3>{0, 0, 1 * UM_IN_ATOMIC_UNITS});

    system_pair_1.diagonalize(diagonalizer);
    system_pair_2.diagonalize(diagonalizer);

    auto eigenenergies_1 = system_pair_1.get_eigenenergies();
    auto eigenenergies_2 = system_pair_2.get_eigenenergies();

    DOCTEST_REQUIRE(eigenenergies_1.size() == eigenenergies_2.size());
    DOCTEST_CHECK(eigenenergies_1.isApprox(eigenenergies_2, 1e-11));
}
} // namespace

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

DOCTEST_TEST_CASE("get amplitudes (via matrix elements) between different pair basis") {
    DiagonalizerEigen<double> diagonalizer;

    Database &database = Database::get_global_instance();
    auto atomic_basis = BasisAtomCreator<double>()
                            .set_species("Rb")
                            .restrict_quantum_number_n(58, 62)
                            .restrict_quantum_number_l(0, 2)
                            .restrict_quantum_number_m(0.5, 0.5)
                            .create(database);

    auto ket = KetAtomCreator()
                   .set_species("Rb")
                   .set_quantum_number_n(60)
                   .set_quantum_number_l(0)
                   .set_quantum_number_m(0.5)
                   .create(database);

    auto perturbed_system1 = SystemAtom<double>(atomic_basis);
    perturbed_system1.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    perturbed_system1.diagonalize(diagonalizer);

    auto perturbed_system2 = SystemAtom<double>(atomic_basis);
    perturbed_system2.set_electric_field({0, 0, 2 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    perturbed_system2.diagonalize(diagonalizer);

    double min_energy = 2 * ket->get_energy() - 20 / HARTREE_IN_GHZ;
    double max_energy = 2 * ket->get_energy() + 20 / HARTREE_IN_GHZ;
    auto perturbed_basis = BasisPairCreator<double>()
                               .add(perturbed_system1)
                               .add(perturbed_system2)
                               .restrict_energy(min_energy, max_energy)
                               .create();
    auto state_in_perturbed_basis = perturbed_basis->get_state(42);

    auto unperturbed_system1 = SystemAtom<double>(atomic_basis);
    auto unperturbed_system2 = SystemAtom<double>(atomic_basis);
    min_energy = 2 * ket->get_energy() - 10 / HARTREE_IN_GHZ;
    max_energy = 2 * ket->get_energy() + 10 / HARTREE_IN_GHZ;
    auto small_unperturbed_basis = BasisPairCreator<double>()
                                       .add(unperturbed_system1)
                                       .add(unperturbed_system2)
                                       .restrict_energy(min_energy, max_energy)
                                       .create();
    DOCTEST_CHECK(perturbed_basis->get_number_of_states() !=
                  small_unperturbed_basis->get_number_of_states());

    auto amplitudes = state_in_perturbed_basis->get_matrix_elements(
        small_unperturbed_basis, OperatorType::IDENTITY, OperatorType::IDENTITY);
    DOCTEST_CHECK(amplitudes.rows() == small_unperturbed_basis->get_number_of_states());
    DOCTEST_CHECK(amplitudes.cols() == state_in_perturbed_basis->get_number_of_states());
}

DOCTEST_TEST_CASE("create a symmetrized BasisPair") {
    auto &database = Database::get_global_instance();
    auto diagonalizer = DiagonalizerEigen<double>();

    auto basis = BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(60, 61)
                     .restrict_quantum_number_l(0, 1)
                     .restrict_quantum_number_m(-0.5, 0.5)
                     .create(database);

    SystemAtom<double> system(basis);
    system.diagonalize(diagonalizer);

    auto canonical_basis = BasisPairCreator<double>().add(system).add(system).create();

    DOCTEST_SUBCASE("restrict permutation parity") {
        auto symmetrized_basis = BasisPairCreator<double>()
                                     .add(system)
                                     .add(system)
                                     .restrict_parity_under_permutation(Parity::ODD)
                                     .create();

        auto expected_basis = canonical_basis->transformed(Transformation<double>(
            build_manual_symmetrizer(canonical_basis, Parity::UNKNOWN, Parity::ODD)));

        check_same_pair_eigenenergies(symmetrized_basis, expected_basis, diagonalizer);
        DOCTEST_CHECK(symmetrized_basis->get_number_of_states() <
                      canonical_basis->get_number_of_states());
    }

    DOCTEST_SUBCASE("restrict inversion parity") {
        auto symmetrized_basis = BasisPairCreator<double>()
                                     .add(system)
                                     .add(system)
                                     .restrict_parity_under_inversion(Parity::ODD)
                                     .create();

        auto expected_basis = canonical_basis->transformed(Transformation<double>(
            build_manual_symmetrizer(canonical_basis, Parity::ODD, Parity::UNKNOWN)));

        check_same_pair_eigenenergies(symmetrized_basis, expected_basis, diagonalizer);
        DOCTEST_CHECK(symmetrized_basis->get_number_of_states() <
                      canonical_basis->get_number_of_states());
    }

    DOCTEST_SUBCASE("restrict permutation parity to EVEN includes identical-state kets") {
        auto symmetrized_basis_even = BasisPairCreator<double>()
                                          .add(system)
                                          .add(system)
                                          .restrict_parity_under_permutation(Parity::EVEN)
                                          .create();

        auto symmetrized_basis_odd = BasisPairCreator<double>()
                                         .add(system)
                                         .add(system)
                                         .restrict_parity_under_permutation(Parity::ODD)
                                         .create();

        // Count kets in the canonical basis where both atoms are in the same state (id1 == id2).
        // Such kets are always permutation-symmetric and must appear in EVEN but not ODD.
        size_t num_diagonal_kets = 0;
        for (const auto &ket : *canonical_basis) {
            auto atomic_states = ket->get_atomic_states();
            DOCTEST_REQUIRE(atomic_states.size() == 2);
            size_t id1 = atomic_states[0]->get_corresponding_ket(0)->get_id_in_database();
            size_t id2 = atomic_states[1]->get_corresponding_ket(0)->get_id_in_database();
            if (id1 == id2) {
                ++num_diagonal_kets;
            }
        }
        DOCTEST_REQUIRE(num_diagonal_kets > 0);

        // EVEN and ODD together partition all canonical kets: off-diagonal pairs contribute one
        // state to each sector, diagonal kets only contribute to EVEN.
        DOCTEST_CHECK(symmetrized_basis_even->get_number_of_states() +
                          symmetrized_basis_odd->get_number_of_states() ==
                      canonical_basis->get_number_of_states());
        DOCTEST_CHECK(symmetrized_basis_odd->get_number_of_states() -
                          symmetrized_basis_even->get_number_of_states() ==
                      num_diagonal_kets);

        // The eigenenergies of the EVEN and ODD bases together must reproduce the
        // eigenenergies of the canonical (non-symmetrized) basis.
        auto system_pair_canonical =
            SystemPair<double>(canonical_basis).set_distance_vector({0, 0, 1 * UM_IN_ATOMIC_UNITS});
        auto system_pair_even = SystemPair<double>(symmetrized_basis_even)
                                    .set_distance_vector({0, 0, 1 * UM_IN_ATOMIC_UNITS});
        auto system_pair_odd = SystemPair<double>(symmetrized_basis_odd)
                                   .set_distance_vector({0, 0, 1 * UM_IN_ATOMIC_UNITS});

        system_pair_canonical.diagonalize(diagonalizer);
        system_pair_even.diagonalize(diagonalizer);
        system_pair_odd.diagonalize(diagonalizer);

        auto canonical_eigenenergies = system_pair_canonical.get_eigenenergies();
        auto even_eigenenergies = system_pair_even.get_eigenenergies();
        auto odd_eigenenergies = system_pair_odd.get_eigenenergies();

        DOCTEST_REQUIRE(canonical_eigenenergies.size() ==
                        even_eigenenergies.size() + odd_eigenenergies.size());

        Eigen::VectorXd combined_eigenenergies(canonical_eigenenergies.size());
        combined_eigenenergies << even_eigenenergies, odd_eigenenergies;
        std::sort(combined_eigenenergies.data(),
                  combined_eigenenergies.data() + combined_eigenenergies.size());
        std::sort(canonical_eigenenergies.data(),
                  canonical_eigenenergies.data() + canonical_eigenenergies.size());

        DOCTEST_CHECK(combined_eigenenergies.isApprox(canonical_eigenenergies, 1e-11));
    }

    DOCTEST_SUBCASE("combine inversion and permutation parity") {
        auto symmetrized_basis = BasisPairCreator<double>()
                                     .add(system)
                                     .add(system)
                                     .restrict_parity_under_inversion(Parity::ODD)
                                     .restrict_parity_under_permutation(Parity::ODD)
                                     .create();

        DOCTEST_CHECK(symmetrized_basis->get_number_of_states() <
                      canonical_basis->get_number_of_states());

        Eigen::SparseMatrix<double, Eigen::ColMajor> coefficients =
            symmetrized_basis->get_coefficients();
        const double inv_sqrt_two = 1 / std::sqrt(2.0);

        for (int state_index = 0; state_index < coefficients.outerSize(); ++state_index) {
            std::vector<std::pair<int, double>> entries;
            for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(coefficients,
                                                                                state_index);
                 it; ++it) {
                auto atomic_states = symmetrized_basis->get_kets()[it.row()]->get_atomic_states();
                DOCTEST_REQUIRE(atomic_states.size() == 2);
                DOCTEST_CHECK(static_cast<int>(atomic_states[0]->get_parity(0)) *
                                  static_cast<int>(atomic_states[1]->get_parity(0)) ==
                              static_cast<int>(Parity::EVEN));
                entries.emplace_back(it.row(), it.value());
            }

            DOCTEST_CHECK(entries.size() >= 1);
            DOCTEST_CHECK(entries.size() <= 2);
            if (entries.size() == 1) {
                DOCTEST_CHECK(entries[0].second == doctest::Approx(1));
            } else {
                DOCTEST_CHECK(std::abs(entries[0].second) == doctest::Approx(inv_sqrt_two));
                DOCTEST_CHECK(std::abs(entries[1].second) == doctest::Approx(inv_sqrt_two));
                DOCTEST_CHECK(entries[0].second == doctest::Approx(entries[1].second));
            }
        }
    }
}

} // namespace pairinteraction
