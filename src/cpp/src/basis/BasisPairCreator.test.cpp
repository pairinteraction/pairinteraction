// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPairCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisAtomCreator.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/enums/SorterType.hpp"
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
#include <limits>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pairinteraction {

constexpr double HARTREE_IN_GHZ = 6579683.920501762;
constexpr double VOLT_PER_CM_IN_ATOMIC_UNITS = 1 / 5.14220675112e9;
constexpr double UM_IN_ATOMIC_UNITS = 1 / 5.29177210544e-5;

namespace {
template <typename Scalar>
std::shared_ptr<const typename BasisAtom<Scalar>::ket_t>
get_dominant_ket(const std::shared_ptr<const BasisAtom<Scalar>> &basis, size_t state_index) {
    Eigen::MatrixX<Scalar> column(
        basis->get_coefficients().col(static_cast<Eigen::Index>(state_index)));
    Eigen::Index ket_index = 0;
    column.cwiseAbs().col(0).maxCoeff(&ket_index);
    return basis->get_ket(static_cast<size_t>(ket_index));
}

// The coefficients of the symmetrized states with respect to the kets of the given basis
template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
build_manually_symmetrized_coefficients(const std::shared_ptr<const BasisPair<Scalar>> &basis,
                                        std::optional<int> parity_under_inversion,
                                        std::optional<int> parity_under_permutation) {
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

            size_t id1 = get_dominant_ket(basis1, idx1)->get_id_in_database();
            size_t id2 = get_dominant_ket(basis2, idx2)->get_id_in_database();

            if (id1 == id2) {
                if (parity_under_inversion == 1 || parity_under_permutation == 1) {
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

            if (id1 < id2) {
                triplets.emplace_back(ket_index, column_index, static_cast<Scalar>(inv_sqrt_two));
            } else {
                int swapped_sign = 0;
                if (parity_under_inversion.has_value() && !parity_under_permutation.has_value()) {
                    swapped_sign = -*parity_under_inversion *
                        static_cast<int>(basis1->get_quantum_number("parity", idx1)) *
                        static_cast<int>(basis2->get_quantum_number("parity", idx2));
                } else {
                    swapped_sign = -*parity_under_permutation;
                }
                triplets.emplace_back(ket_index, column_index,
                                      static_cast<Scalar>(swapped_sign * inv_sqrt_two));
            }
        }
    }

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients(
        static_cast<Eigen::Index>(basis->get_number_of_kets()), state_index);
    coefficients.setFromTriplets(triplets.begin(), triplets.end());
    return coefficients;
}

template <typename Scalar>
std::shared_ptr<const BasisPair<Scalar>>
build_pair_basis(std::shared_ptr<const BasisAtom<Scalar>> basis1,
                 std::shared_ptr<const BasisAtom<Scalar>> basis2) {
    return BasisPairCreator<Scalar>(std::move(basis1), std::move(basis2)).create();
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
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 58, 62)
                     .restrict_quantum_number("l", 0, 2)
                     .create();
    SystemAtom<double> system(basis);
    system.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});

    DiagonalizerEigen<double> diagonalizer;
    system.diagonalize(diagonalizer);

    // Get energy window for a two-atom basis
    auto ket = KetAtomCreator(database)
                   .set_species("Rb")
                   .set_quantum_number("n", 60)
                   .set_quantum_number("l", 0)
                   .set_quantum_number("m", 0.5)
                   .create();
    double min_energy = 2 * ket->get_energy() - 3 / HARTREE_IN_GHZ;
    double max_energy = 2 * ket->get_energy() + 3 / HARTREE_IN_GHZ;

    // Create two-atom bases
    auto basis_pair_a =
        pairinteraction::BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
            .restrict_energy(min_energy, max_energy)
            .restrict_quantum_number("m", 1, 1)
            .create();
    auto basis_pair_b =
        pairinteraction::BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
            .restrict_energy(min_energy, max_energy)
            .restrict_quantum_number("m", 1, 1)
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
        DOCTEST_CHECK(*ket1a == *ket1b);
        DOCTEST_CHECK(*ket2a == *ket2b);
    }

    DOCTEST_SUBCASE("check overlap") {
        auto basis_ket = BasisAtomCreator<double>(database).add_ket(ket).create();
        auto basis_pair_ket = build_pair_basis<double>(basis_ket, basis_ket);
        Eigen::RowVectorXd amplitudes =
            basis_pair_a
                ->get_matrix_elements(basis_pair_ket, OperatorType::IDENTITY,
                                      OperatorType::IDENTITY, 0, 0)
                .row(0);
        auto overlaps = amplitudes.cwiseAbs2().eval();

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

DOCTEST_TEST_CASE("merge rejects differently transformed atomic bases") {
    auto &database = Database::get_global_instance();
    auto atomic_basis = BasisAtomCreator<double>(database)
                            .set_species("Rb")
                            .restrict_quantum_number("n", 60, 60)
                            .restrict_quantum_number("l", 0, 1)
                            .restrict_quantum_number("m", 0.5, 0.5)
                            .create();

    auto system_a = SystemAtom<double>(atomic_basis);
    system_a.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    system_a.diagonalize(DiagonalizerEigen<double>());

    auto system_b = SystemAtom<double>(atomic_basis);
    system_b.set_electric_field({0, 0, 2 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    system_b.diagonalize(DiagonalizerEigen<double>());

    auto pair_basis_a =
        BasisPairCreator<double>(system_a.get_eigenbasis(), system_a.get_eigenbasis()).create();
    auto pair_basis_b =
        BasisPairCreator<double>(system_b.get_eigenbasis(), system_b.get_eigenbasis()).create();

    DOCTEST_CHECK_THROWS_WITH_AS(pair_basis_a->merge(pair_basis_b),
                                 doctest::Contains("Cannot merge two pair bases"),
                                 std::invalid_argument);
}

DOCTEST_TEST_CASE("get matrix elements in the pair basis") {
    DiagonalizerEigen<double> diagonalizer;

    // Create single-atom system
    Database &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 58, 62)
                     .restrict_quantum_number("l", 0, 2)
                     .create();
    SystemAtom<double> system(basis);
    system.set_electric_field({0, 0, 10 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    system.diagonalize(diagonalizer);

    // Get energy window for a two-atom basis
    auto ket = KetAtomCreator(database)
                   .set_species("Rb")
                   .set_quantum_number("n", 60)
                   .set_quantum_number("l", 0)
                   .set_quantum_number("m", 0.5)
                   .create();
    double min_energy = 2 * ket->get_energy() - 3 / HARTREE_IN_GHZ;
    double max_energy = 2 * ket->get_energy() + 3 / HARTREE_IN_GHZ;

    // Create two-atom system
    auto basis_pair_unperturbed =
        pairinteraction::BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
            .restrict_energy(min_energy, max_energy)
            .restrict_quantum_number("m", 1, 1)
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
        auto atomic_states = basis_pair_unperturbed->get_kets()[0]->get_atomic_states();
        auto basis_pair_ket_pair = build_pair_basis<double>(atomic_states[0], atomic_states[1]);
        Eigen::RowVectorXd matrix_elements_ket_pair =
            basis_pair_unperturbed
                ->get_matrix_elements(basis_pair_ket_pair, OperatorType::ELECTRIC_DIPOLE,
                                      OperatorType::ELECTRIC_DIPOLE, 0, 0)
                .row(0);
        DOCTEST_CHECK(matrix_elements_ket_pair.size() ==
                      basis_pair_unperturbed->get_number_of_states());

        {
            // The basis is symmetrized, so the first ket is a superposition of the states. Because
            // the symmetrization is a real orthogonal transformation, the coefficients of this
            // superposition are given by the first row of the coefficient matrix.
            Eigen::RowVectorXd ref =
                basis_pair_unperturbed->get_coefficients().row(0) * matrix_elements_all;
            DOCTEST_CHECK(ref.isApprox(matrix_elements_ket_pair, 1e-11));
        }

        // <basis x basis|d0d0|basis_pair>
        auto basis_pair_product = build_pair_basis<double>(basis, basis);
        auto matrix_elements_product = basis_pair->get_matrix_elements(
            basis_pair_product, OperatorType::ELECTRIC_DIPOLE, OperatorType::ELECTRIC_DIPOLE, 0, 0);
        DOCTEST_CHECK(matrix_elements_product.rows() ==
                      basis->get_number_of_states() * basis->get_number_of_states());
        DOCTEST_CHECK(matrix_elements_product.cols() == basis_pair->get_number_of_states());

        // <ket,ket|d0d0|basis_pair>
        auto basis_ket = BasisAtomCreator<double>(database).add_ket(ket).create();
        auto basis_pair_ket = build_pair_basis<double>(basis_ket, basis_ket);
        auto matrix_elements_ket = basis_pair->get_matrix_elements(
            basis_pair_ket, OperatorType::ELECTRIC_DIPOLE, OperatorType::ELECTRIC_DIPOLE, 0, 0);
        DOCTEST_CHECK(matrix_elements_ket.rows() == 1);
        DOCTEST_CHECK(matrix_elements_ket.cols() == basis_pair->get_number_of_states());
    }

    DOCTEST_SUBCASE("check matrix elements") {
        // energy
        auto hamiltonian = basis_pair->get_matrix_elements(basis_pair, OperatorType::ENERGY,
                                                           OperatorType::IDENTITY, 0, 0);
        hamiltonian += basis_pair->get_matrix_elements(basis_pair, OperatorType::IDENTITY,
                                                       OperatorType::ENERGY, 0, 0);

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
    auto atomic_basis = BasisAtomCreator<double>(database)
                            .set_species("Rb")
                            .restrict_quantum_number("n", 58, 62)
                            .restrict_quantum_number("l", 0, 2)
                            .restrict_quantum_number("m", 0.5, 0.5)
                            .create();

    auto ket = KetAtomCreator(database)
                   .set_species("Rb")
                   .set_quantum_number("n", 60)
                   .set_quantum_number("l", 0)
                   .set_quantum_number("m", 0.5)
                   .create();

    auto perturbed_system1 = SystemAtom<double>(atomic_basis);
    perturbed_system1.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    perturbed_system1.diagonalize(diagonalizer);

    auto perturbed_system2 = SystemAtom<double>(atomic_basis);
    perturbed_system2.set_electric_field({0, 0, 2 * VOLT_PER_CM_IN_ATOMIC_UNITS});
    perturbed_system2.diagonalize(diagonalizer);

    double min_energy = 2 * ket->get_energy() - 20 / HARTREE_IN_GHZ;
    double max_energy = 2 * ket->get_energy() + 20 / HARTREE_IN_GHZ;
    auto perturbed_basis = BasisPairCreator<double>(perturbed_system1.get_eigenbasis(),
                                                    perturbed_system2.get_eigenbasis())
                               .restrict_energy(min_energy, max_energy)
                               .create();
    auto state_in_perturbed_basis = perturbed_basis->get_state(42);

    min_energy = 2 * ket->get_energy() - 10 / HARTREE_IN_GHZ;
    max_energy = 2 * ket->get_energy() + 10 / HARTREE_IN_GHZ;
    auto small_unperturbed_basis = BasisPairCreator<double>(atomic_basis, atomic_basis)
                                       .restrict_energy(min_energy, max_energy)
                                       .create();
    DOCTEST_CHECK(perturbed_basis->get_number_of_states() !=
                  small_unperturbed_basis->get_number_of_states());

    auto amplitudes = state_in_perturbed_basis->get_matrix_elements(
        small_unperturbed_basis, OperatorType::IDENTITY, OperatorType::IDENTITY, 0, 0);
    DOCTEST_CHECK(amplitudes.rows() == small_unperturbed_basis->get_number_of_states());
    DOCTEST_CHECK(amplitudes.cols() == state_in_perturbed_basis->get_number_of_states());
}

namespace {
// The atomic basis and its eigenbasis that are shared by the tests of the symmetrized pair basis.
// As both atoms share the same BasisAtom, the pair basis can be symmetrized.
struct SymmetrizedBasisPairFixture {
    DiagonalizerEigen<double> diagonalizer;
    std::shared_ptr<const BasisAtom<double>> basis =
        BasisAtomCreator<double>(Database::get_global_instance())
            .set_species("Rb")
            .restrict_quantum_number("n", 60, 61)
            .restrict_quantum_number("l", 0, 1)
            .restrict_quantum_number("m", -0.5, 0.5)
            .create();
    SystemAtom<double> system{basis};
    std::shared_ptr<const BasisPair<double>> canonical_basis;

    SymmetrizedBasisPairFixture() {
        system.diagonalize(diagonalizer);

        // As both atoms share the same BasisAtom, the created basis is symmetrized. We undo the
        // symmetrization to compare against a manual symmetrization.
        canonical_basis = BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                              .create()
                              ->canonicalized();
    }
};
} // namespace

DOCTEST_TEST_CASE_FIXTURE(SymmetrizedBasisPairFixture, "create a symmetrized BasisPair") {
    DOCTEST_SUBCASE("restrict permutation parity") {
        auto symmetrized_basis =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                .restrict_quantum_number("parity_under_permutation", -1, -1)
                .create();

        auto symmetrized_coefficients =
            build_manually_symmetrized_coefficients(canonical_basis, std::nullopt, -1);
        auto expected_basis = canonical_basis->get_modified_copy(
            symmetrized_coefficients,
            std::vector<double>(static_cast<size_t>(symmetrized_coefficients.cols()),
                                std::numeric_limits<double>::max()),
            {});

        check_same_pair_eigenenergies(symmetrized_basis, expected_basis, diagonalizer);
        DOCTEST_CHECK(symmetrized_basis->get_number_of_states() <
                      canonical_basis->get_number_of_states());

        // A state of the symmetrized basis has the energy of the kets it comprises, which share the
        // same energy
        DOCTEST_REQUIRE(symmetrized_basis->has_energy());
        const auto &coefficients = symmetrized_basis->get_coefficients();
        for (Eigen::Index row = 0; row < coefficients.outerSize(); ++row) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(coefficients, row);
                 it; ++it) {
                DOCTEST_CHECK(symmetrized_basis->get_energy(static_cast<size_t>(it.col())) ==
                              symmetrized_basis->get_ket(static_cast<size_t>(row))->get_energy());
            }
        }
    }

    DOCTEST_SUBCASE("restrict inversion parity") {
        auto symmetrized_basis =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                .restrict_quantum_number("parity_under_inversion", -1, -1)
                .create();

        auto symmetrized_coefficients =
            build_manually_symmetrized_coefficients(canonical_basis, -1, std::nullopt);
        auto expected_basis = canonical_basis->get_modified_copy(
            symmetrized_coefficients,
            std::vector<double>(static_cast<size_t>(symmetrized_coefficients.cols()),
                                std::numeric_limits<double>::max()),
            {});

        check_same_pair_eigenenergies(symmetrized_basis, expected_basis, diagonalizer);
        DOCTEST_CHECK(symmetrized_basis->get_number_of_states() <
                      canonical_basis->get_number_of_states());
    }

    DOCTEST_SUBCASE("restrict permutation parity to EVEN includes identical-state kets") {
        auto symmetrized_basis_even =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                .restrict_quantum_number("parity_under_permutation", 1, 1)
                .create();

        auto symmetrized_basis_odd =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                .restrict_quantum_number("parity_under_permutation", -1, -1)
                .create();

        // Count kets in the canonical basis where both atoms are in the same state (id1 == id2).
        // Such kets are always permutation-symmetric and must appear in EVEN but not ODD.
        size_t num_diagonal_kets = 0;
        for (const auto &ket : *canonical_basis) {
            auto atomic_states = ket->get_atomic_states();
            DOCTEST_REQUIRE(atomic_states.size() == 2);
            size_t id1 = get_dominant_ket(atomic_states[0], 0)->get_id_in_database();
            size_t id2 = get_dominant_ket(atomic_states[1], 0)->get_id_in_database();
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

    DOCTEST_SUBCASE("restrict the product of parities") {
        for (bool symmetrization_enabled : {false, true}) {
            auto restricted_basis =
                BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                    .set_symmetrization_enabled(symmetrization_enabled)
                    .restrict_quantum_number("product_of_parities", -1, -1)
                    .create();
            DOCTEST_REQUIRE(restricted_basis->get_number_of_states() > 0);
            for (const auto &ket : restricted_basis->get_kets()) {
                DOCTEST_CHECK(ket->get_quantum_number("product_of_parities") == -1);
            }
        }
    }

    DOCTEST_SUBCASE("the symmetrization is enabled by default if possible") {
        // A second, independently constructed system represents a different atom. Even though it
        // is built from the same basis, it is a distinct object, so symmetrization is impossible.
        SystemAtom<double> system_other(basis);
        system_other.diagonalize(diagonalizer);

        // For the same BasisAtom, the pair basis is symmetrized by default ...
        auto symmetrized_basis =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis()).create();
        DOCTEST_CHECK_FALSE(symmetrized_basis->is_canonical());
        DOCTEST_CHECK(symmetrized_basis->has_quantum_number("parity_under_permutation"));

        // ... unless the symmetrization is disabled
        auto product_basis =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                .set_symmetrization_enabled(false)
                .create();
        DOCTEST_CHECK(product_basis->is_canonical());
        DOCTEST_CHECK_FALSE(product_basis->has_quantum_number("parity_under_permutation"));

        // For different BasisAtom objects, the pair basis cannot be symmetrized
        auto different_basis =
            BasisPairCreator<double>(system.get_eigenbasis(), system_other.get_eigenbasis())
                .create();
        DOCTEST_CHECK(different_basis->is_canonical());
        DOCTEST_CHECK_THROWS_AS(
            BasisPairCreator<double>(system.get_eigenbasis(), system_other.get_eigenbasis())
                .set_symmetrization_enabled(true),
            std::invalid_argument);

        // The parities are only well-defined for a symmetrized pair basis, so that they cannot be
        // restricted, not even by a range comprising both parities
        DOCTEST_CHECK_THROWS_AS(
            BasisPairCreator<double>(system.get_eigenbasis(), system_other.get_eigenbasis())
                .restrict_quantum_number("parity_under_permutation", -1, -1)
                .create(),
            std::invalid_argument);
        DOCTEST_CHECK_THROWS_AS(
            BasisPairCreator<double>(system.get_eigenbasis(), system_other.get_eigenbasis())
                .restrict_quantum_number("parity_under_inversion", -1, 1)
                .create(),
            std::invalid_argument);
    }
}

DOCTEST_TEST_CASE_FIXTURE(SymmetrizedBasisPairFixture,
                          "BasisPairCreator validates the quantum number restrictions") {
    DOCTEST_SUBCASE("quantum number restrictions are validated") {
        BasisPairCreator<double> creator(system.get_eigenbasis(), system.get_eigenbasis());

        // Parities must be +1 or -1
        DOCTEST_CHECK_THROWS_AS(creator.restrict_quantum_number("parity_under_inversion", 0, 0),
                                std::invalid_argument);
        DOCTEST_CHECK_THROWS_AS(creator.restrict_quantum_number("parity_under_permutation", 2, 2),
                                std::invalid_argument);

        // The quantum number m must be an integer or half-integer
        DOCTEST_CHECK_THROWS_AS(creator.restrict_quantum_number("m", 0.25, 1),
                                std::invalid_argument);

        // The minimum must not exceed the maximum
        DOCTEST_CHECK_THROWS_AS(creator.restrict_quantum_number("parity_under_inversion", 1, -1),
                                std::invalid_argument);

        // Only quantum numbers of pair states can be restricted
        DOCTEST_CHECK_THROWS_AS(creator.restrict_quantum_number("n", 60, 60),
                                std::invalid_argument);
        DOCTEST_CHECK_THROWS_AS(creator.restrict_quantum_number("parity", 1, 1),
                                std::invalid_argument);
    }

    DOCTEST_SUBCASE("a range comprising both parities does not restrict the parity") {
        auto unrestricted_basis =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis()).create();
        auto basis_with_both_parities =
            BasisPairCreator<double>(system.get_eigenbasis(), system.get_eigenbasis())
                .restrict_quantum_number("parity_under_inversion", -1, 1)
                .restrict_quantum_number("parity_under_permutation", -1, 1)
                .create();
        DOCTEST_CHECK(basis_with_both_parities->get_coefficients().isApprox(
            unrestricted_basis->get_coefficients()));
    }
}

DOCTEST_TEST_CASE_FIXTURE(SymmetrizedBasisPairFixture,
                          "label the states of the symmetrized pair basis by their parities") {
    DOCTEST_SUBCASE("label every state by its parities") {
        auto eigenbasis = system.get_eigenbasis();

        // Following https://doi.org/10.1088/1361-6455/aa743a, a state with the label p is an
        // eigenstate of the permutation and inversion operators with the eigenvalue -p
        size_t number_of_states = 0;
        for (int parity_under_permutation : {1, -1}) {
            for (int parity_under_inversion : {1, -1}) {
                auto restricted_basis =
                    BasisPairCreator<double>(eigenbasis, eigenbasis)
                        .restrict_quantum_number("parity_under_permutation",
                                                 parity_under_permutation, parity_under_permutation)
                        .restrict_quantum_number("parity_under_inversion", parity_under_inversion,
                                                 parity_under_inversion)
                        .create();
                number_of_states += restricted_basis->get_number_of_states();

                const auto number_of_kets =
                    static_cast<Eigen::Index>(restricted_basis->get_number_of_kets());
                Eigen::SparseMatrix<double, Eigen::RowMajor> permutation(number_of_kets,
                                                                         number_of_kets);
                Eigen::SparseMatrix<double, Eigen::RowMajor> inversion(number_of_kets,
                                                                       number_of_kets);
                for (size_t idx1 = 0; idx1 < eigenbasis->get_number_of_states(); ++idx1) {
                    for (size_t idx2 = 0; idx2 < eigenbasis->get_number_of_states(); ++idx2) {
                        int col = restricted_basis->get_ket_index_from_tuple(idx1, idx2);
                        if (col < 0) {
                            continue;
                        }
                        int row = restricted_basis->get_ket_index_from_tuple(idx2, idx1);
                        DOCTEST_REQUIRE(row >= 0);
                        permutation.insert(row, col) = 1;
                        inversion.insert(row, col) =
                            eigenbasis->get_quantum_number("parity", idx1) *
                            eigenbasis->get_quantum_number("parity", idx2);
                    }
                }

                Eigen::MatrixXd coefficients = restricted_basis->get_coefficients();
                for (size_t i = 0; i < restricted_basis->get_number_of_states(); ++i) {
                    DOCTEST_CHECK(restricted_basis->get_quantum_number(
                                      "parity_under_permutation", i) == parity_under_permutation);
                    DOCTEST_CHECK(restricted_basis->get_quantum_number(
                                      "parity_under_inversion", i) == parity_under_inversion);
                    Eigen::VectorXd state = coefficients.col(static_cast<Eigen::Index>(i));
                    DOCTEST_CHECK(state.norm() == doctest::Approx(1));
                    DOCTEST_CHECK(
                        (permutation * state).isApprox(-parity_under_permutation * state, 1e-12));
                    DOCTEST_CHECK(
                        (inversion * state).isApprox(-parity_under_inversion * state, 1e-12));
                }
            }
        }

        // The sectors partition the basis, and canonicalizing drops the labels
        DOCTEST_CHECK(number_of_states == canonical_basis->get_number_of_states());
        DOCTEST_CHECK_FALSE(canonical_basis->has_quantum_number("parity_under_permutation"));
        DOCTEST_CHECK_FALSE(canonical_basis->has_quantum_number("parity_under_inversion"));
    }

    DOCTEST_SUBCASE("only label the parity under permutation if the atomic parity is undefined") {
        SystemAtom<double> system_in_field(basis);
        system_in_field.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
        system_in_field.diagonalize(diagonalizer);
        auto eigenbasis = system_in_field.get_eigenbasis();

        auto symmetrized_basis = BasisPairCreator<double>(eigenbasis, eigenbasis).create();
        DOCTEST_CHECK(symmetrized_basis->has_quantum_number("parity_under_permutation"));
        DOCTEST_CHECK_FALSE(symmetrized_basis->has_quantum_number("parity_under_inversion"));

        DOCTEST_CHECK_THROWS_AS(BasisPairCreator<double>(eigenbasis, eigenbasis)
                                    .restrict_quantum_number("parity_under_inversion", -1, -1)
                                    .create(),
                                std::invalid_argument);
    }
}

DOCTEST_TEST_CASE("BasisPairCreator checks the energies of the BasisAtom") {
    auto &database = Database::get_global_instance();
    auto basis = BasisAtomCreator<double>(database)
                     .set_species("Rb")
                     .restrict_quantum_number("n", 60, 61)
                     .restrict_quantum_number("l", 0, 1)
                     .create();
    double min_energy = 2 * basis->get_energy(0);
    double max_energy = 2 * basis->get_energy(basis->get_number_of_states() - 1);

    DOCTEST_SUBCASE("a canonical basis sorted by energy can be used") {
        DOCTEST_CHECK_NOTHROW(BasisPairCreator<double>(basis, basis).create());
        DOCTEST_CHECK_NOTHROW(BasisPairCreator<double>(basis, basis)
                                  .restrict_energy(min_energy, max_energy)
                                  .create());
    }

    DOCTEST_SUBCASE("the states must have well-defined energies for an energy restriction") {
        auto basis_without_energies = basis->get_modified_copy(
            basis->get_coefficients(),
            std::vector<double>(basis->get_number_of_states(), std::numeric_limits<double>::max()),
            basis->get_quantum_numbers_of_states());
        DOCTEST_CHECK_THROWS_WITH_AS(BasisPairCreator<double>(basis, basis_without_energies)
                                         .restrict_energy(min_energy, max_energy)
                                         .create(),
                                     doctest::Contains("well-defined energies"),
                                     std::invalid_argument);

        // Without an energy restriction, the pair states do not have a well-defined energy
        auto pair_basis = BasisPairCreator<double>(basis, basis_without_energies).create();
        DOCTEST_CHECK(pair_basis->get_number_of_states() ==
                      basis->get_number_of_states() * basis->get_number_of_states());
        DOCTEST_CHECK_FALSE(pair_basis->has_energy());
    }

    DOCTEST_SUBCASE("the states must be sorted by energy for an energy restriction") {
        auto basis_sorted_by_m = basis->sorted(basis->get_sorter({SorterType::QUANTUM_NUMBER_M}));
        DOCTEST_CHECK_THROWS_WITH_AS(BasisPairCreator<double>(basis_sorted_by_m, basis)
                                         .restrict_energy(min_energy, max_energy)
                                         .create(),
                                     doctest::Contains("sorted by energy"), std::invalid_argument);
        DOCTEST_CHECK_NOTHROW(BasisPairCreator<double>(basis_sorted_by_m, basis).create());
    }

    DOCTEST_SUBCASE("the eigenbasis of a diagonalized system can be used") {
        SystemAtom<double> system(basis);
        system.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
        system.diagonalize(DiagonalizerEigen<double>());
        auto eigenbasis = system.get_eigenbasis();
        auto pair_basis = BasisPairCreator<double>(eigenbasis, eigenbasis).create();

        // The energies of the pair states are the sums of the eigenenergies of the atoms
        for (const auto &ket : *pair_basis) {
            auto atomic_states = ket->get_atomic_states();
            DOCTEST_CHECK(
                ket->get_energy() ==
                doctest::Approx(atomic_states[0]->get_energy(0) + atomic_states[1]->get_energy(0)));
        }

        // The eigenbasis of a system that was not sorted by energy cannot be used with an energy
        // restriction
        SystemAtom<double> unsorted_system(basis);
        unsorted_system.set_electric_field({0, 0, 1 * VOLT_PER_CM_IN_ATOMIC_UNITS});
        unsorted_system.diagonalize(DiagonalizerEigen<double>(), {}, {}, 1e-6, false);
        DOCTEST_CHECK_THROWS_WITH_AS(
            BasisPairCreator<double>(unsorted_system.get_eigenbasis(), eigenbasis)
                .restrict_energy(min_energy, max_energy)
                .create(),
            doctest::Contains("sorted by energy"), std::invalid_argument);
    }
}

} // namespace pairinteraction
