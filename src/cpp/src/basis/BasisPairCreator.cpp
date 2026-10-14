// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPairCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/utils/TaskControl.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
BasisPairCreator<Scalar>::BasisPairCreator(std::shared_ptr<const BasisAtom<Scalar>> basis1,
                                           std::shared_ptr<const BasisAtom<Scalar>> basis2)
    : basis1(std::move(basis1)), basis2(std::move(basis2)),
      symmetrization_enabled(this->basis1 == this->basis2) {}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::restrict_energy(real_t min, real_t max) {
    range_energy = {min, max};
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &
BasisPairCreator<Scalar>::restrict_quantum_number(const std::string &name, real_t min, real_t max) {
    if (name == "m") {
        if (2 * min != std::rint(2 * min) || 2 * max != std::rint(2 * max)) {
            throw std::invalid_argument("Quantum number m must be an integer or half-integer.");
        }
    } else if (name == "product_of_parities" || name == "parity_under_inversion" ||
               name == "parity_under_permutation") {
        if ((min != 1 && min != -1) || (max != 1 && max != -1)) {
            throw std::invalid_argument("The parity must be +1 or -1.");
        }
    } else {
        throw std::invalid_argument("The quantum number " + name + " cannot be restricted.");
    }
    quantum_number_ranges[name] = {min, max};
    return *this;
}

template <typename Scalar>
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::set_symmetrization_enabled(bool value) {
    // The pair states can only be symmetrized if both atoms share the same BasisAtom, so that a
    // one-atom state can be identified across both atoms by its state index
    if (value && basis1 != basis2) {
        throw std::invalid_argument(
            "The symmetrization requires the same BasisAtom for both atoms, because it is only "
            "defined for two identical atoms.");
    }
    symmetrization_enabled = value;
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisPair<Scalar>> BasisPairCreator<Scalar>::create() const {
    set_task_status("Constructing pair basis...");

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();
    static constexpr real_t unknown = std::numeric_limits<real_t>::max();

    // Sorting is required for the binary search of the energetically allowed range
    if (range_energy.is_finite() &&
        (!basis1->is_sorted_by_energy() || !basis2->is_sorted_by_energy())) {
        throw std::invalid_argument(
            "The energy can only be restricted if the states of the BasisAtom have well-defined "
            "energies and are sorted by energy. Consider using the eigenbasis of a diagonalized "
            "SystemAtom.");
    }

    // ---------------------------------------------------------------------------------------------
    // Determine which quantum numbers of the pair states are well-defined
    // ---------------------------------------------------------------------------------------------

    const bool has_quantum_number_m =
        basis1->has_quantum_number("m") && basis2->has_quantum_number("m");
    const bool has_product_of_parities =
        basis1->has_quantum_number("parity") && basis2->has_quantum_number("parity");

    const std::unordered_map<std::string, bool> is_well_defined{
        {"m", has_quantum_number_m},
        {"product_of_parities", has_product_of_parities},
        {"parity_under_permutation", symmetrization_enabled},
        {"parity_under_inversion", symmetrization_enabled && has_product_of_parities},
    };
    for (const auto &[name, range] : quantum_number_ranges) {
        if (!is_well_defined.at(name)) {
            throw std::invalid_argument(
                "The quantum number " + name +
                " must not be restricted because it is not well-defined. Note that the parities "
                "under inversion and permutation are only well-defined if the symmetrization is "
                "enabled. The parity under inversion and the product of parities additionally "
                "require the one-atom states to have a well-defined parity.");
        }
    }

    // ---------------------------------------------------------------------------------------------
    // Define helpers for adding kets and states
    // ---------------------------------------------------------------------------------------------

    using quantum_numbers_of_state_t = std::unordered_map<std::string, double>;

    ketvec_t kets;
    kets.reserve(basis1->get_number_of_states() * basis2->get_number_of_states());
    typename basis_t::map_range_t state_index1_to_state_index_range2;
    state_index1_to_state_index_range2.reserve(basis1->get_number_of_states());
    typename basis_t::map_indices_t state_indices_to_ket_index;

    auto add_ket = [&](size_t idx1, size_t idx2, real_t energy,
                       const quantum_numbers_of_state_t &quantum_numbers) {
        auto ket_index = static_cast<Eigen::Index>(kets.size());
        state_indices_to_ket_index.try_emplace(std::vector<size_t>{idx1, idx2}, kets.size());
        kets.emplace_back(std::make_shared<ket_t>(
            typename ket_t::Private(), std::initializer_list<size_t>{idx1, idx2},
            std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>>{basis1, basis2}, energy,
            quantum_numbers));
        return ket_index;
    };

    std::vector<Eigen::Triplet<Scalar>> triplets;
    std::vector<real_t> energy_of_states;
    typename basis_t::quantum_numbers_t quantum_numbers_of_states;
    Eigen::Index num_states = 0;

    // Get the quantum numbers of the symmetrized state that shares the quantum numbers of its
    // kets and has the given parity under permutation. Following
    // https://doi.org/10.1088/1361-6455/aa743a, its parity under inversion is the product of the
    // parity under permutation and the product of parities.
    auto get_quantum_numbers_of_state = [&](quantum_numbers_of_state_t quantum_numbers,
                                            int parity_under_permutation) {
        quantum_numbers["parity_under_permutation"] = parity_under_permutation;
        if (has_product_of_parities) {
            quantum_numbers["parity_under_inversion"] =
                parity_under_permutation * quantum_numbers.at("product_of_parities");
        }
        return quantum_numbers;
    };

    // Check whether all quantum numbers of a state are within their ranges. Quantum numbers
    // without a restriction are always allowed.
    auto is_allowed_state = [&](const quantum_numbers_of_state_t &quantum_numbers) {
        for (const auto &[name, value] : quantum_numbers) {
            auto it = quantum_number_ranges.find(name);
            if (it == quantum_number_ranges.end()) {
                continue;
            }
            if (value < it->second.min() - numerical_precision ||
                value > it->second.max() + numerical_precision) {
                return false;
            }
        }
        return true;
    };

    // Add the state that is the superposition of the given kets (row index and coefficient). It
    // shares the energy of the kets.
    auto add_state = [&](std::initializer_list<std::pair<Eigen::Index, Scalar>> components,
                         real_t energy, const quantum_numbers_of_state_t &quantum_numbers) {
        for (const auto &[row_index, coefficient] : components) {
            triplets.emplace_back(row_index, num_states, coefficient);
        }
        energy_of_states.push_back(energy);
        for (const auto &[name, value] : quantum_numbers) {
            quantum_numbers_of_states[name].push_back(static_cast<real_t>(value));
        }
        ++num_states;
    };

    // ---------------------------------------------------------------------------------------------
    // Construct the kets and the states
    // ---------------------------------------------------------------------------------------------

    // States without a well-defined energy are marked by the maximum value
    auto get_energies = [](const auto &basis) {
        std::vector<real_t> energies(basis->get_number_of_states(), unknown);
        if (basis->has_energy()) {
            for (size_t state_index = 0; state_index < energies.size(); ++state_index) {
                energies[state_index] = basis->get_energy(state_index);
            }
        }
        return energies;
    };
    const std::vector<real_t> eigenenergies1 = get_energies(basis1);
    const std::vector<real_t> eigenenergies2 = get_energies(basis2);
    const real_t *eigenenergies2_begin = eigenenergies2.data();
    const real_t *eigenenergies2_end = eigenenergies2_begin + eigenenergies2.size();

    const auto inverse_sqrt_two = static_cast<Scalar>(1 / std::sqrt(2.0));

    for (size_t idx1 = 0; idx1 < eigenenergies1.size(); ++idx1) {
        set_task_status("Constructing pair basis...");

        // Get the range of the second index for which the pair energy is allowed.
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

        // If symmetrizing, the kets |a, b> and |b, a> are handled together for a <= b
        if (symmetrization_enabled) {
            min = std::max(min, idx1);
        }

        for (size_t idx2 = min; idx2 < max; ++idx2) {
            const real_t energy =
                (eigenenergies1[idx1] == unknown || eigenenergies2[idx2] == unknown)
                ? unknown
                : eigenenergies1[idx1] + eigenenergies2[idx2];

            // Get the quantum numbers of the ket. They are symmetric under the exchange of the
            // atoms, so that all checks below are symmetric, too.
            quantum_numbers_of_state_t quantum_numbers;
            if (has_quantum_number_m) {
                quantum_numbers["m"] =
                    basis1->get_quantum_number("m", idx1) + basis2->get_quantum_number("m", idx2);
            }
            if (has_product_of_parities) {
                quantum_numbers["product_of_parities"] =
                    basis1->get_quantum_number("parity", idx1) *
                    basis2->get_quantum_number("parity", idx2);
            }

            // Without symmetrization, the ket itself is a state of the basis
            if (!symmetrization_enabled) {
                if (is_allowed_state(quantum_numbers)) {
                    auto row_index = add_ket(idx1, idx2, energy, quantum_numbers);
                    add_state({{row_index, 1}}, energy, quantum_numbers);
                }
                continue;
            }

            // The ket |a, a> is symmetric under the exchange of the atoms. Following
            // https://doi.org/10.1088/1361-6455/aa743a, it has the parity under permutation -1.
            if (idx1 == idx2) {
                auto quantum_numbers_of_state = get_quantum_numbers_of_state(quantum_numbers, -1);
                if (is_allowed_state(quantum_numbers_of_state)) {
                    auto row_index = add_ket(idx1, idx2, energy, quantum_numbers);
                    add_state({{row_index, 1}}, energy, quantum_numbers_of_state);
                }
                continue;
            }

            // Following https://doi.org/10.1088/1361-6455/aa743a, the symmetrized states are
            // (|a, b> - f |b, a>) / sqrt(2) for a < b with the parity under permutation f. Here,
            // it is a = idx1 and b = idx2. The kets are only added if one of the states is allowed.
            const auto quantum_numbers_plus = get_quantum_numbers_of_state(quantum_numbers, 1);
            const auto quantum_numbers_minus = get_quantum_numbers_of_state(quantum_numbers, -1);
            const bool is_allowed_plus = is_allowed_state(quantum_numbers_plus);
            const bool is_allowed_minus = is_allowed_state(quantum_numbers_minus);
            if (!is_allowed_plus && !is_allowed_minus) {
                continue;
            }
            auto row_index_ab = add_ket(idx1, idx2, energy, quantum_numbers);
            auto row_index_ba = add_ket(idx2, idx1, energy, quantum_numbers);
            if (is_allowed_plus) {
                add_state({{row_index_ba, inverse_sqrt_two}, {row_index_ab, -inverse_sqrt_two}},
                          energy, quantum_numbers_plus);
            }
            if (is_allowed_minus) {
                add_state({{row_index_ba, inverse_sqrt_two}, {row_index_ab, inverse_sqrt_two}},
                          energy, quantum_numbers_minus);
            }
        }
    }
    kets.shrink_to_fit();

    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> coefficients(
        static_cast<Eigen::Index>(kets.size()), num_states);
    coefficients.setFromTriplets(triplets.begin(), triplets.end());

    return std::make_shared<const basis_t>(
        typename basis_t::Private(), std::move(kets), std::move(state_index1_to_state_index_range2),
        std::move(state_indices_to_ket_index), basis1, basis2, std::move(coefficients),
        std::move(energy_of_states), std::move(quantum_numbers_of_states));
}

// Explicit instantiations
template class BasisPairCreator<double>;
template class BasisPairCreator<std::complex<double>>;
} // namespace pairinteraction
