#include "pairinteraction/basis/BasisCombinedCreator.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisCombined.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/enums/TransformationType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetCombined.hpp"
#include "pairinteraction/system/SystemAtom.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fmt/core.h>
#include <limits>
#include <memory>

namespace pairinteraction {
template <typename Scalar>
BasisCombinedCreator<Scalar> &
BasisCombinedCreator<Scalar>::add(const SystemAtom<Scalar> &system_atom) {
    if (!system_atom.is_diagonal()) {
        throw std::invalid_argument("The system must be diagonalized before it can be added.");
    }
    systems_atom.push_back(system_atom);
    return *this;
}

template <typename Scalar>
BasisCombinedCreator<Scalar> &BasisCombinedCreator<Scalar>::restrict_energy(real_t min,
                                                                            real_t max) {
    range_energy = {min, max};
    return *this;
}

template <typename Scalar>
BasisCombinedCreator<Scalar> &BasisCombinedCreator<Scalar>::restrict_quantum_number_m(real_t min,
                                                                                      real_t max) {
    range_quantum_number_m = {min, max};
    return *this;
}

template <typename Scalar>
std::shared_ptr<const BasisCombined<Scalar>> BasisCombinedCreator<Scalar>::create() const {
    if (systems_atom.size() != 2) {
        throw std::invalid_argument(
            "Two SystemAtom must be added before creating the combined basis.");
    }

    real_t numerical_precision = 10 * std::numeric_limits<real_t>::epsilon();

    // Sort the states, which are eigenstates, by their energy
    auto system1 = systems_atom[0].get();
    auto system2 = systems_atom[1].get();
    system1.transform(system1.get_sorter({TransformationType::SORT_BY_ENERGY}));
    system2.transform(system2.get_sorter({TransformationType::SORT_BY_ENERGY}));

    // Construct the canonical basis that contains all combined states with allowed energies and
    // quantum numbers
    auto basis1 = system1.get_basis();
    auto basis2 = system2.get_basis();
    auto eigenvalues1 = system1.get_eigenvalues();
    auto eigenvalues2 = system2.get_eigenvalues();
    real_t *eigenvalues2_begin = eigenvalues2.data();
    real_t *eigenvalues2_end = eigenvalues2_begin + eigenvalues2.size();

    ketvec_t kets;
    kets.reserve(eigenvalues1.size() * eigenvalues2.size());

    typename basis_t::map_range_t map_range_of_state_index2;
    map_range_of_state_index2.reserve(eigenvalues1.size());

    typename basis_t::map_indices_t state_indices_to_ket_index;

    // Loop only over states with an allowed energy
    size_t ket_index = 0;
    for (size_t idx1 = 0; idx1 < static_cast<size_t>(eigenvalues1.size()); ++idx1) {
        // Get the energetically allowed range of the second index
        size_t min = 0;
        size_t max = eigenvalues2.size();
        if (range_energy.is_finite()) {
            real_t min_val2 = range_energy.min() - eigenvalues1[idx1];
            real_t max_val2 = range_energy.max() - eigenvalues1[idx1];
            min = std::distance(eigenvalues2_begin,
                                std::lower_bound(eigenvalues2_begin, eigenvalues2_end, min_val2));
            max = std::distance(eigenvalues2_begin,
                                std::upper_bound(eigenvalues2_begin, eigenvalues2_end, max_val2));
        }
        map_range_of_state_index2.emplace(idx1, typename basis_t::range_t(min, max));

        // Loop over the energetically allowed range of the second index
        for (size_t idx2 = min; idx2 < max; ++idx2) {
            // Get energy
            const real_t energy = eigenvalues1[idx1] + eigenvalues2[idx2];
            assert(!range_energy.is_finite() ||
                   (energy >= range_energy.min() && energy <= range_energy.max()));

            // Get labels for the kets with largest overlap
            auto label1 = basis1->get_corresponding_ket(idx1)->get_label();
            auto label2 = basis2->get_corresponding_ket(idx2)->get_label();

            // Create a combined state
            auto ket = std::make_shared<ket_t>(
                typename ket_t::Private(), std::initializer_list<size_t>{idx1, idx2},
                std::initializer_list<std::string>{label1, label2},
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

            // Store the combined state as a ket
            kets.emplace_back(std::move(ket));

            // Store the ket index of the combined state
            state_indices_to_ket_index.emplace(std::vector<size_t>{idx1, idx2}, ket_index++);
        }
    }

    kets.shrink_to_fit();

    auto now = std::chrono::high_resolution_clock::now();
    std::uint64_t nanoseconds =
        std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
    std::string id_of_kets = fmt::format("{:016x}", nanoseconds);

    return std::make_shared<basis_t>(typename basis_t::Private(), std::move(kets),
                                     std::move(id_of_kets), std::move(map_range_of_state_index2),
                                     std::move(state_indices_to_ket_index), basis1, basis2);
}

// Explicit instantiations
template class BasisCombinedCreator<float>;
template class BasisCombinedCreator<double>;
template class BasisCombinedCreator<std::complex<float>>;
template class BasisCombinedCreator<std::complex<double>>;
} // namespace pairinteraction
