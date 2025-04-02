// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/BasisPairCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <algorithm>
#include <limits>
#include <memory>

namespace pairinteraction {
template <typename Scalar>
BasisPairCreator<Scalar>::BasisPairCreator() : product_of_parities(Parity::UNKNOWN) {}

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
BasisPairCreator<Scalar> &BasisPairCreator<Scalar>::restrict_product_of_parities(Parity value) {
    product_of_parities = value;
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisPair<Scalar>> BasisPairCreator<Scalar>::create() const {
    if (systems_atom.size() != 2) {
        throw std::invalid_argument("Two SystemAtom must be added before creating the BasisPair.");
    }

    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

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

    // Loop only over states with an allowed energy
    size_t ket_index = 0;
    for (size_t idx1 = 0; idx1 < static_cast<size_t>(eigenenergies1.size()); ++idx1) {
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

            // Check the parity of the sum of the parities
            if (product_of_parities != Parity::UNKNOWN) {
                if (static_cast<int>(basis1->get_parity(idx1)) *
                        static_cast<int>(basis2->get_parity(idx2)) !=
                    static_cast<int>(product_of_parities)) {
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

            // Store the ket index
            state_indices_to_ket_index.emplace(std::vector<size_t>{idx1, idx2}, ket_index++);
        }
    }

    kets.shrink_to_fit();

    return std::make_shared<basis_t>(typename basis_t::Private(), std::move(kets),
                                     std::move(map_range_of_state_index2),
                                     std::move(state_indices_to_ket_index), basis1, basis2);
}

// Explicit instantiations
template class BasisPairCreator<double>;
template class BasisPairCreator<std::complex<double>>;
} // namespace pairinteraction
