// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPairCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"
#include "pairinteraction/utils/TaskControl.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::add(const SystemAtom<Scalar> &system_atom) {
    // The system must be diagonalized and its eigenstates sorted by energy.
    // Sorting is required for the binary search of the energetically allowed range in create().
    // By default, System::diagonalize ensures this.
    if (!system_atom.is_diagonal_and_sorted_by_energy()) {
        throw std::invalid_argument(
            "The system must be diagonalized and sorted by energy before it can be added. "
            "Consider calling diagonalize() on the SystemAtom which also sorts the eigenstates.");
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
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::restrict_parity_under_inversion(int value) {
    if (value != 1 && value != -1) {
        throw std::invalid_argument("The parity must be +1 or -1.");
    }
    parity_under_inversion = value;
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::restrict_parity_under_permutation(int value) {
    if (value != 1 && value != -1) {
        throw std::invalid_argument("The parity must be +1 or -1.");
    }
    parity_under_permutation = value;
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::set_symmetrization_enabled(bool enable) {
    symmetrization_enabled = enable;
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisPair<Scalar>> BasisPairCreator<Scalar>::create() const {
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    set_task_status("Constructing pair basis...");

    // ---------------------------------------------------------------------------------------------
    // Check the systems
    // ---------------------------------------------------------------------------------------------

    if (systems_atom.size() != 2) {
        throw std::invalid_argument("Two SystemAtom must be added before creating the BasisPair.");
    }

    // Only references to the systems are stored, so a system might have been changed since add()
    for (const auto &system_atom : systems_atom) {
        if (!system_atom.get().is_diagonal_and_sorted_by_energy()) {
            throw std::invalid_argument(
                "The systems must still be diagonalized and sorted by energy when the BasisPair is "
                "created. Do not change a SystemAtom after it has been added.");
        }
    }

    const auto &system1 = systems_atom[0].get();
    const auto &system2 = systems_atom[1].get();
    auto basis1 = system1.get_basis();
    auto basis2 = system2.get_basis();
    auto eigenenergies1 = system1.get_eigenenergies();
    auto eigenenergies2 = system2.get_eigenenergies();

    // ---------------------------------------------------------------------------------------------
    // Check which restrictions and quantum numbers are well-defined
    // ---------------------------------------------------------------------------------------------

    // Symmetrization is only defined for two identical atoms. Requiring the same SystemAtom to be
    // added twice ensures that a one-atom state can be identified across both atoms by its state
    // index. If it has not been explicitly enabled or disabled, symmetrization is applied whenever
    // this is possible.
    const bool is_symmetrization_possible = &system1 == &system2;
    const bool has_parity_restriction =
        parity_under_inversion.has_value() || parity_under_permutation.has_value();
    if (symmetrization_enabled == false && has_parity_restriction) {
        throw std::invalid_argument(
            "Parity restrictions require symmetrization, which has been disabled.");
    }
    if (!is_symmetrization_possible && (symmetrization_enabled == true || has_parity_restriction)) {
        throw std::invalid_argument(
            "Symmetrization and parity restrictions require the same SystemAtom to be added twice, "
            "because symmetrization is only defined for two identical atoms.");
    }
    const bool do_symmetrization = symmetrization_enabled.value_or(is_symmetrization_possible);

    // The quantum number m of a pair state is only well-defined if it is well-defined for both
    // atoms
    const bool has_quantum_number_m =
        basis1->has_quantum_number("m") && basis2->has_quantum_number("m");
    if (!has_quantum_number_m && range_quantum_number_m.is_finite()) {
        throw std::invalid_argument(
            "The quantum number m must not be restricted because it is not well-defined.");
    }

    // The parity under inversion of a symmetrized pair state is only well-defined if the parities
    // of the one-atom states are well-defined
    const bool has_parity_under_inversion =
        do_symmetrization && basis1->has_quantum_number("parity");
    if (!has_parity_under_inversion && parity_under_inversion.has_value()) {
        throw std::invalid_argument(
            "The parity under inversion must not be restricted because it is not well-defined, "
            "as the one-atom states do not have a well-defined parity.");
    }

    auto get_product_of_parities = [&](size_t idx1, size_t idx2) {
        return static_cast<int>(basis1->get_quantum_number("parity", idx1)) *
            static_cast<int>(basis2->get_quantum_number("parity", idx2));
    };

    // ---------------------------------------------------------------------------------------------
    // Construct the kets, i.e. all product states with allowed energies and quantum numbers
    // ---------------------------------------------------------------------------------------------

    // Kets that cannot contribute to any symmetrized pair state of the requested parities are
    // skipped right away:
    // - Following https://doi.org/10.1088/1361-6455/aa743a, pair states |a, a> are always of odd
    //   parity, so they do not contribute if an even parity is requested.
    // - The parity under inversion is the product of the parity under permutation and the parities
    //   of the one-atom states. Thus, if both are restricted, only kets with a certain product of
    //   parities contribute.
    // Both criteria are symmetric under the exchange of the atoms.
    const bool exclude_identical_states =
        do_symmetrization && (parity_under_inversion == 1 || parity_under_permutation == 1);
    std::optional<int> required_product_of_parities;
    if (parity_under_inversion.has_value() && parity_under_permutation.has_value()) {
        required_product_of_parities = *parity_under_inversion * *parity_under_permutation;
    }
    auto can_contribute = [&](size_t idx1, size_t idx2) {
        if (exclude_identical_states && idx1 == idx2) {
            return false;
        }
        return !required_product_of_parities.has_value() ||
            get_product_of_parities(idx1, idx2) == *required_product_of_parities;
    };

    ketvec_t kets;
    kets.reserve(eigenenergies1.size() * eigenenergies2.size());
    typename basis_t::map_range_t state_index1_to_state_index_range2;
    state_index1_to_state_index_range2.reserve(eigenenergies1.size());
    typename basis_t::map_indices_t state_indices_to_ket_index;

    const real_t *eigenenergies2_begin = eigenenergies2.data();
    const real_t *eigenenergies2_end = eigenenergies2_begin + eigenenergies2.size();
    for (size_t idx1 = 0; idx1 < static_cast<size_t>(eigenenergies1.size()); ++idx1) {
        set_task_status("Constructing pair basis kets...");

        // Get the range of indices of the second atom's states for which the pair energy is
        // allowed. As the pair energy itself is compared to the energy window, |a, b> is allowed
        // if and only if |b, a> is allowed.
        size_t min = 0;
        size_t max = eigenenergies2.size();
        if (range_energy.is_finite()) {
            const real_t energy1 = eigenenergies1[idx1];
            min = std::distance(
                eigenenergies2_begin,
                std::partition_point(eigenenergies2_begin, eigenenergies2_end, [&](real_t energy2) {
                    return energy1 + energy2 < range_energy.min();
                }));
            max = std::distance(
                eigenenergies2_begin,
                std::partition_point(eigenenergies2_begin, eigenenergies2_end, [&](real_t energy2) {
                    return energy1 + energy2 <= range_energy.max();
                }));
        }
        state_index1_to_state_index_range2.try_emplace(idx1, typename basis_t::range_t(min, max));

        for (size_t idx2 = min; idx2 < max; ++idx2) {
            if (!can_contribute(idx1, idx2)) {
                continue;
            }

            std::unordered_map<std::string, double> quantum_numbers;
            if (has_quantum_number_m) {
                const real_t m =
                    basis1->get_quantum_number("m", idx1) + basis2->get_quantum_number("m", idx2);
                if (range_quantum_number_m.is_finite() &&
                    (m < range_quantum_number_m.min() - numerical_precision ||
                     m > range_quantum_number_m.max() + numerical_precision)) {
                    continue;
                }
                quantum_numbers["m"] = m;
            }

            const real_t energy = eigenenergies1[idx1] + eigenenergies2[idx2];
            assert(!range_energy.is_finite() ||
                   (energy >= range_energy.min() && energy <= range_energy.max()));

            state_indices_to_ket_index.try_emplace(std::vector<size_t>{idx1, idx2}, kets.size());
            kets.emplace_back(std::make_shared<ket_t>(
                typename ket_t::Private(), std::initializer_list<size_t>{idx1, idx2},
                std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>>{basis1, basis2},
                energy, std::move(quantum_numbers)));
        }
    }
    kets.shrink_to_fit();

    if (!do_symmetrization) {
        return std::make_shared<basis_t>(typename basis_t::Private(), std::move(kets),
                                         std::move(state_index1_to_state_index_range2),
                                         std::move(state_indices_to_ket_index), basis1, basis2);
    }

    // ---------------------------------------------------------------------------------------------
    // Construct the symmetrized pair states
    // ---------------------------------------------------------------------------------------------

    // The symmetrized states are
    //   |a, a> with parity_under_permutation = parity_under_inversion = -1, and
    //   (|a, b> - p |b, a>) / sqrt(2) for a > b with parity_under_permutation = p and
    //   parity_under_inversion = p * parity(a) * parity(b).
    // Every state is labeled by both parities, but only states matching the restrictions are kept.
    set_task_status("Symmetrizing pair basis...");

    std::vector<Eigen::Triplet<Scalar>> coefficient_triplets;
    coefficient_triplets.reserve(kets.size());
    typename basis_t::quantum_numbers_of_states_t quantum_numbers_of_states;
    for (const auto &[label, name] : basis_t::sorter_type_to_quantum_number_name) {
        quantum_numbers_of_states[name].reserve(kets.size());
    }
    Eigen::Index number_of_states = 0;

    // Add a state that is the superposition of the given kets (row index and coefficient) if it
    // has the requested parities. The ket passed by reference provides all other quantum numbers,
    // which are shared by all kets contributing to the state.
    auto add_state = [&](const ket_t &ket, int permutation,
                         std::initializer_list<std::pair<Eigen::Index, Scalar>> components) {
        std::optional<int> inversion;
        if (has_parity_under_inversion) {
            inversion =
                permutation * get_product_of_parities(ket.atomic_indices[0], ket.atomic_indices[1]);
        }
        if ((parity_under_permutation.has_value() && permutation != *parity_under_permutation) ||
            (parity_under_inversion.has_value() && inversion != *parity_under_inversion)) {
            return;
        }

        for (const auto &[row_index, coefficient] : components) {
            coefficient_triplets.emplace_back(row_index, number_of_states, coefficient);
        }

        for (const auto &[label, name] : basis_t::sorter_type_to_quantum_number_name) {
            real_t quantum_number = std::numeric_limits<real_t>::max();
            if (label == SorterType::PARITY_UNDER_PERMUTATION) {
                quantum_number = permutation;
            } else if (label == SorterType::PARITY_UNDER_INVERSION) {
                if (inversion.has_value()) {
                    quantum_number = *inversion;
                }
            } else if (ket.has_quantum_number(name)) {
                quantum_number = static_cast<real_t>(ket.get_quantum_number(name));
            }
            quantum_numbers_of_states[name].push_back(quantum_number);
        }

        ++number_of_states;
    };

    const auto inverse_sqrt_two = static_cast<Scalar>(1 / std::sqrt(2.0));
    for (Eigen::Index row_index = 0; row_index < static_cast<Eigen::Index>(kets.size());
         ++row_index) {
        const auto &ket = *kets[row_index];
        const size_t idx1 = ket.atomic_indices[0];
        const size_t idx2 = ket.atomic_indices[1];

        if (idx1 == idx2) {
            add_state(ket, -1, {{row_index, 1}});
            continue;
        }

        // The kets |a, b> and |b, a> are combined when encountering the ket with a > b
        if (idx1 < idx2) {
            continue;
        }

        // The partner ket exists because all filters above are symmetric under the exchange of
        // the atoms
        auto it = state_indices_to_ket_index.find(std::vector<size_t>{idx2, idx1});
        if (it == state_indices_to_ket_index.end()) {
            throw std::logic_error("The partner ket required for symmetrization is missing.");
        }
        const auto partner_row_index = static_cast<Eigen::Index>(it->second);

        for (int permutation : {1, -1}) {
            add_state(ket, permutation,
                      {{row_index, inverse_sqrt_two},
                       {partner_row_index, static_cast<Scalar>(-permutation) * inverse_sqrt_two}});
        }
    }

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients(
        static_cast<Eigen::Index>(kets.size()), number_of_states);
    coefficients.setFromTriplets(coefficient_triplets.begin(), coefficient_triplets.end());

    return std::make_shared<basis_t>(typename basis_t::Private(), std::move(kets),
                                     std::move(coefficients), std::move(quantum_numbers_of_states),
                                     std::move(state_index1_to_state_index_range2),
                                     std::move(state_indices_to_ket_index), basis1, basis2);
}

// Explicit instantiations
template class BasisPairCreator<double>;
template class BasisPairCreator<std::complex<double>>;
} // namespace pairinteraction
