// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPairCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/TaskControl.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace pairinteraction {
template <typename Scalar>
BasisPairCreator<Scalar>::BasisPairCreator()
    : parity_under_inversion(Parity::UNKNOWN), parity_under_permutation(Parity::UNKNOWN) {}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::add(const SystemAtom<Scalar> &system_atom) {
    if (!system_atom.is_diagonal()) {
        throw std::invalid_argument("The system must be diagonalized before it can be added.");
    }
    systems_atom.push_back(system_atom);
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::restrict_energy(real_t min, real_t max) {
    range_energy = {min, max};
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::restrict_quantum_number_m(real_t min,
                                                                              real_t max) {
    range_quantum_number_m = {min, max};
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::restrict_parity_under_inversion(Parity value) {
    parity_under_inversion = value;
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &
BasisPairCreator<Scalar>::restrict_parity_under_permutation(Parity value) {
    parity_under_permutation = value;
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisPair<Scalar>> BasisPairCreator<Scalar>::create() const {
    set_task_status("Constructing pair basis...");

    if (systems_atom.size() != 2) {
        throw std::invalid_argument("Two SystemAtom must be added before creating the BasisPair.");
    }

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();
    const bool has_symmetry_restriction =
        parity_under_inversion != Parity::UNKNOWN || parity_under_permutation != Parity::UNKNOWN;

    Parity inferred_product_of_parities = Parity::UNKNOWN;
    if (parity_under_inversion != Parity::UNKNOWN && parity_under_permutation != Parity::UNKNOWN) {
        inferred_product_of_parities = static_cast<Parity>(
            static_cast<int>(parity_under_inversion) * static_cast<int>(parity_under_permutation));
    }

    // Sort the states, which are eigenstates, by their energy
    auto system1 = systems_atom[0].get();
    auto system2 = systems_atom[1].get();
    system1.transform(system1.get_sorter({TransformationType::SORT_BY_ENERGY}));
    system2.transform(system2.get_sorter({TransformationType::SORT_BY_ENERGY}));

    // Construct the canonical basis that contains all KetPair objects with allowed energies and
    // quantum numbers
    auto basis1 = system1.get_basis();
    auto basis2 = system2.get_basis();
    auto eigenenergies1 = system1.get_eigenenergies();
    auto eigenenergies2 = system2.get_eigenenergies();
    real_t *eigenenergies2_begin = eigenenergies2.data();
    real_t *eigenenergies2_end = eigenenergies2_begin + eigenenergies2.size();

    ketvec_t kets;
    kets.reserve(eigenenergies1.size() * eigenenergies2.size());

    typename basis_t::map_range_t map_range_of_state_index2;
    map_range_of_state_index2.reserve(eigenenergies1.size());

    typename basis_t::map_indices_t state_indices_to_ket_index;

    Eigen::Index state_index = 0;
    std::unordered_map<std::array<size_t, 2>, Eigen::Index, utils::hash<std::array<size_t, 2>>>
        ket_ids2state_index;
    std::vector<Eigen::Triplet<Scalar>> transformation_triplets;

    if (has_symmetry_restriction) {
        transformation_triplets.reserve(eigenenergies1.size() * eigenenergies2.size());
    }

    // Loop only over states with an allowed energy
    size_t ket_index = 0;
    for (size_t idx1 = 0; idx1 < static_cast<size_t>(eigenenergies1.size()); ++idx1) {
        set_task_status("Constructing pair basis...");

        // Get the energetically allowed range of the second index
        size_t min = 0;
        size_t max = eigenenergies2.size();
        if (range_energy.is_finite()) {
            real_t min_val2 = range_energy.min() - eigenenergies1[idx1];
            real_t max_val2 = range_energy.max() - eigenenergies1[idx1];
            min =
                std::distance(eigenenergies2_begin,
                              std::lower_bound(eigenenergies2_begin, eigenenergies2_end, min_val2));
            max =
                std::distance(eigenenergies2_begin,
                              std::upper_bound(eigenenergies2_begin, eigenenergies2_end, max_val2));
        }
        map_range_of_state_index2.emplace(idx1, typename basis_t::range_t(min, max));

        // Loop over the energetically allowed range of the second index
        for (size_t idx2 = min; idx2 < max; ++idx2) {
            // Get energy
            const real_t energy = eigenenergies1[idx1] + eigenenergies2[idx2];
            assert(!range_energy.is_finite() ||
                   (energy >= range_energy.min() && energy <= range_energy.max()));

            // Check the parity of the product of the parities
            if (inferred_product_of_parities != Parity::UNKNOWN) {
                if (static_cast<int>(basis1->get_parity(idx1)) *
                        static_cast<int>(basis2->get_parity(idx2)) !=
                    static_cast<int>(inferred_product_of_parities)) {
                    continue;
                }
            }

            // Create a KetPair object
            auto ket = std::make_shared<ket_t>(
                typename ket_t::Private(), std::initializer_list<size_t>{idx1, idx2},
                std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>>{basis1, basis2},
                energy);

            // Check the quantum number m
            if (ket->has_quantum_number_m()) {
                if (range_quantum_number_m.is_finite() &&
                    (ket->get_quantum_number_m() <
                         range_quantum_number_m.min() - numerical_precision ||
                     ket->get_quantum_number_m() >
                         range_quantum_number_m.max() + numerical_precision)) {
                    continue;
                }
            } else if (range_quantum_number_m.is_finite()) {
                throw std::invalid_argument(
                    "The quantum number m must not be restricted because it is not well-defined.");
            }

            // Store the KetPair object as a ket
            kets.emplace_back(std::move(ket));

            // Store the ket index after all symmetry contributions have been recorded
            state_indices_to_ket_index.emplace(std::vector<size_t>{idx1, idx2}, ket_index);

            // Store a symmetrized coefficient matrix only if inversion or permutation symmetry is
            // specified
            auto row_index = static_cast<Eigen::Index>(ket_index++);
            if (!has_symmetry_restriction) {
                continue;
            }

            // TODO: on the long run, the states should have IDs and not only the kets. This is in
            // particular of importance if get_corresponding_ket is not working because fields lead
            // to states with contributions from many kets.
            size_t id1 = basis1->get_corresponding_ket(idx1)->get_id_in_database();
            size_t id2 = basis2->get_corresponding_ket(idx2)->get_id_in_database();

            if (id1 == id2 && parity_under_inversion != Parity::EVEN &&
                parity_under_permutation != Parity::EVEN) {
                // Obtain the column index of the symmetrized state
                Eigen::Index column_index = state_index++;

                // A pair state with identical one-atom states contributes with coefficient one
                // if allowed by both symmetries.
                transformation_triplets.emplace_back(row_index, column_index, 1);
            } else if (id1 > id2) {
                // Obtain the column index of the symmetrized state
                std::array<size_t, 2> ordered_ids{std::max(id1, id2), std::min(id1, id2)};
                auto [column_iterator, inserted] =
                    ket_ids2state_index.try_emplace(ordered_ids, state_index);
                if (inserted) {
                    ++state_index;
                }
                Eigen::Index column_index = column_iterator->second;

                // We choose pair states with id1 > id2 to contribute with the coefficient 1/sqrt(2)
                // and put the phase in the contribution of the pair states with id1 < id2.
                transformation_triplets.emplace_back(row_index, column_index, 1 / std::sqrt(2.0));
            } else if (id1 < id2) {
                // Obtain the column index of the symmetrized state
                std::array<size_t, 2> ordered_ids{std::max(id1, id2), std::min(id1, id2)};
                auto [column_iterator, inserted] =
                    ket_ids2state_index.try_emplace(ordered_ids, state_index);
                if (inserted) {
                    ++state_index;
                }
                Eigen::Index column_index = column_iterator->second;

                // Obtain the phase of the contribution to the symmetrized state
                int phase{0};
                if (parity_under_permutation != Parity::UNKNOWN) {
                    // If both inversion and permutation are restricted, the earlier filter on the
                    // product of parities already guarantees that the phases in the inversion and
                    // permutation symmetric states are the same.
                    phase = -static_cast<int>(parity_under_permutation);
                } else {
                    // If only inversion is restricted, the phase of the contribution to the
                    // symmetrized state is determined by the product of the parities of the
                    // one-atom states and the specified inversion parity.
                    phase = -static_cast<int>(parity_under_inversion) *
                        static_cast<int>(basis1->get_parity(idx1)) *
                        static_cast<int>(basis2->get_parity(idx2));
                }

                transformation_triplets.emplace_back(row_index, column_index,
                                                     phase * 1 / std::sqrt(2.0));
            }
        }
    }

    kets.shrink_to_fit();

    std::shared_ptr<const basis_t> basis = std::make_shared<basis_t>(
        typename basis_t::Private(), std::move(kets), std::move(map_range_of_state_index2),
        std::move(state_indices_to_ket_index), basis1, basis2);

    if (!has_symmetry_restriction) {
        return basis;
    }

    transformation_triplets.shrink_to_fit();

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> transformation_matrix(
        basis->get_number_of_states(), state_index);
    transformation_matrix.setFromTriplets(transformation_triplets.begin(),
                                          transformation_triplets.end());

    Eigen::Matrix<real_t, Eigen::Dynamic, 1> sum_of_squared_coefficients =
        transformation_matrix.cwiseAbs2().transpose() *
        Eigen::Matrix<real_t, Eigen::Dynamic, 1>::Ones(transformation_matrix.rows());
    for (Eigen::Index column_index = 0; column_index < sum_of_squared_coefficients.size();
         ++column_index) {
        if (std::abs(sum_of_squared_coefficients[column_index] - 1) > numerical_precision) {
            throw std::invalid_argument(
                "The basis could not be symmetrized. This likely means that the specified parity "
                "restrictions are invalid for the given one-atom systems.");
        }
    }

    // TODO: on the long run, construct the coefficient matrix directly
    auto transformation = Transformation<Scalar>(std::move(transformation_matrix));
    return basis->transformed(transformation);
}

// Explicit instantiations
template class BasisPairCreator<double>;
template class BasisPairCreator<std::complex<double>>;
} // namespace pairinteraction
