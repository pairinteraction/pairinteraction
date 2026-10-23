// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/basis/Basis.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/basis/BasisPair.hpp"
#include "pairinteraction/enums/SorterType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/ket/KetPair.hpp"
#include "pairinteraction/utils/TaskControl.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace pairinteraction {

template <typename Scalar>
class BasisAtom;

namespace {
// The states can only be labeled by quantum numbers by which they can be sorted
template <typename Derived>
void check_quantum_number_name(const std::string &name) {
    if (std::none_of(Derived::sorter_type_to_quantum_number_name.begin(),
                     Derived::sorter_type_to_quantum_number_name.end(),
                     [&name](const auto &entry) { return entry.second == name; })) {
        throw std::invalid_argument("The states cannot be labeled by the quantum number " + name +
                                    ".");
    }
}
} // namespace

template <typename Derived>
void Basis<Derived>::perform_blocks_checks(const std::set<SorterType> &unique_labels) const {
    // Currently, since states in a basis dont store the energy, they cannot be sorted by energy.
    if (unique_labels.contains(SorterType::ENERGY)) {
        throw std::invalid_argument("Blocks cannot be obtained by the energy. Note that sorting "
                                    "a system by the energy is supported nevertheless.");
    }

    // Check if the states are labeled by the requested labels
    std::vector<const std::vector<real_t> *> quantum_numbers_of_labels;
    quantum_numbers_of_labels.reserve(unique_labels.size());
    for (const auto &label : unique_labels) {
        quantum_numbers_of_labels.push_back(&get_quantum_numbers_of_states(label));
    }

    // Check if the states are sorted by the requested labels, i.e., if states that share the same
    // labels are contiguous. Otherwise, states belonging to the same block would be scattered over
    // several blocks and couplings between them would be lost.
    // The quantum numbers are (half-)integers, thus doubling them yields an exact integer
    using label_t = std::vector<int64_t>;
    auto label_of = [&](Eigen::Index i) {
        label_t label;
        label.reserve(quantum_numbers_of_labels.size());
        for (const auto *quantum_numbers : quantum_numbers_of_labels) {
            label.push_back(std::llround(2 * (*quantum_numbers)[i]));
        }
        return label;
    };

    // Collect the labels of the contiguous blocks of states
    std::vector<label_t> labels_of_blocks;
    for (Eigen::Index i = 0; i < coefficients.cols(); ++i) {
        label_t label = label_of(i);
        if (labels_of_blocks.empty() || labels_of_blocks.back() != label) {
            labels_of_blocks.push_back(label);
        }
    }

    // If the states are sorted, no two blocks share the same labels
    std::sort(labels_of_blocks.begin(), labels_of_blocks.end());
    if (std::adjacent_find(labels_of_blocks.begin(), labels_of_blocks.end()) !=
        labels_of_blocks.end()) {
        throw std::invalid_argument("The states are not sorted by the requested labels.");
    }
}

template <typename Derived>
Basis<Derived>::Basis(ketvec_t &&kets,
                      Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &&coefficients,
                      energies_t &&energy_of_states, quantum_numbers_t &&quantum_numbers_of_states)
    : kets(std::move(kets)), coefficients(std::move(coefficients)),
      quantum_numbers_of_states(std::move(quantum_numbers_of_states)),
      energy_of_states(std::move(energy_of_states)) {
    if (this->kets.empty()) {
        throw std::invalid_argument("The basis must contain at least one element.");
    }
    check_states();
}

template <typename Derived>
void Basis<Derived>::check_states() const {
    auto num_states = static_cast<size_t>(coefficients.cols());

    // Neither infinite values nor the maximum value are accepted as energies or quantum numbers,
    // so that no value can be mistaken for marking an energy or quantum number as not well-defined
    auto is_valid_value = [](real_t value) {
        return std::isfinite(value) && std::abs(value) < std::numeric_limits<real_t>::max();
    };

    if (static_cast<size_t>(coefficients.rows()) != kets.size()) {
        throw std::invalid_argument(
            "The number of rows of the coefficient matrix must match the number of kets.");
    }

    // The energy is either given for all states or not well-defined for the basis
    if (energy_of_states.has_value()) {
        if (energy_of_states->size() != num_states) {
            throw std::invalid_argument("The number of energies must match the number of states.");
        }
        if (!std::ranges::all_of(*energy_of_states, is_valid_value)) {
            throw std::invalid_argument(
                "The energies of the states must be finite and smaller than the maximum value.");
        }
    }

    // A quantum number is either given for all states or not well-defined for the basis
    for (const auto &[name, quantum_numbers] : quantum_numbers_of_states) {
        check_quantum_number_name<Derived>(name);
        if (quantum_numbers.size() != num_states) {
            throw std::invalid_argument("The number of values of the quantum number " + name +
                                        " must match the number of states.");
        }

        // The quantum numbers are (half-)integers, so that doubling them yields the exact integers
        // by which the states are sorted and block-diagonalized
        for (real_t quantum_number : quantum_numbers) {
            real_t doubled = 2 * quantum_number;
            if (!is_valid_value(quantum_number) || !std::isfinite(doubled) ||
                doubled != std::round(doubled)) {
                throw std::invalid_argument("The quantum number " + name +
                                            " must be an integer or a half-integer.");
            }
        }
    }
}

template <typename Derived>
typename Basis<Derived>::StatesData Basis<Derived>::get_canonical_states(const ketvec_t &kets) {
    StatesData states;
    auto n = static_cast<Eigen::Index>(kets.size());

    // Make the i-th state equal to the i-th ket
    states.coefficients.resize(n, n);
    states.coefficients.setIdentity();

    // Assign each state the quantum numbers of the ket it is equal to. A quantum number is only
    // assigned if all kets have it.
    for (const auto &entry : Derived::sorter_type_to_quantum_number_name) {
        const std::string &name = entry.second;
        if (!std::ranges::all_of(
                kets, [&name](const auto &ket) { return ket->has_quantum_number(name); })) {
            continue;
        }
        auto &quantum_numbers = states.quantum_numbers_of_states[name];
        quantum_numbers.reserve(kets.size());
        for (const auto &ket : kets) {
            quantum_numbers.push_back(static_cast<real_t>(ket->get_quantum_number(name)));
        }
    }

    // Assign each state the energy of the ket it is equal to. The energy is only assigned if all
    // kets have it.
    if (std::ranges::all_of(kets, [](const auto &ket) { return ket->has_energy(); })) {
        auto &energies = states.energy_of_states.emplace();
        energies.reserve(kets.size());
        for (const auto &ket : kets) {
            energies.push_back(static_cast<real_t>(ket->get_energy()));
        }
    }

    return states;
}

template <typename Derived>
bool Basis<Derived>::has_quantum_number(const std::string &name) const {
    return quantum_numbers_of_states.contains(name);
}

template <typename Derived>
bool Basis<Derived>::has_energy() const {
    return energy_of_states.has_value();
}

template <typename Derived>
bool Basis<Derived>::is_sorted_by_energy() const {
    return has_energy() && std::ranges::is_sorted(*energy_of_states);
}

template <typename Derived>
const Derived &Basis<Derived>::derived() const {
    return static_cast<const Derived &>(*this);
}

template <typename Derived>
const std::vector<typename Basis<Derived>::real_t> &
Basis<Derived>::get_quantum_numbers_of_states(SorterType label) const {
    auto it = Derived::sorter_type_to_quantum_number_name.find(label);
    if (it == Derived::sorter_type_to_quantum_number_name.end()) {
        throw std::invalid_argument(
            "The states of this basis cannot be labeled by the requested sorter type.");
    }
    const std::string &name = it->second;
    if (!has_quantum_number(name)) {
        throw std::invalid_argument(
            "States cannot be labeled and thus not sorted by the quantum number " + name + ".");
    }
    return quantum_numbers_of_states.at(name);
}

template <typename Derived>
const typename Basis<Derived>::ketvec_t &Basis<Derived>::get_kets() const {
    return kets;
}

template <typename Derived>
const Eigen::SparseMatrix<typename Basis<Derived>::scalar_t, Eigen::RowMajor> &
Basis<Derived>::get_coefficients() const {
    return coefficients;
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_quantum_number(const std::string &name,
                                                                   size_t state_index) const {
    // A quantum number by which the states are not labeled is not well-defined
    auto it = quantum_numbers_of_states.find(name);
    if (it == quantum_numbers_of_states.end()) {
        check_quantum_number_name<Derived>(name);
        throw std::invalid_argument("The state does not have a well-defined quantum number " +
                                    name + ".");
    }
    return it->second.at(state_index);
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_energy(size_t state_index) const {
    if (!has_energy()) {
        throw std::invalid_argument("The state does not have a well-defined energy.");
    }
    return energy_of_states->at(state_index);
}

template <typename Derived>
const typename Basis<Derived>::quantum_numbers_t &
Basis<Derived>::get_quantum_numbers_of_states() const {
    return quantum_numbers_of_states;
}

template <typename Derived>
std::shared_ptr<const Derived> Basis<Derived>::get_state(size_t state_index) const {
    // Create a basis that consists of the single requested state, without copying the other states
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> coefficients_of_state =
        coefficients.col(static_cast<Eigen::Index>(state_index));
    quantum_numbers_t quantum_numbers_of_state;
    for (const auto &[name, quantum_numbers] : quantum_numbers_of_states) {
        quantum_numbers_of_state[name] = {quantum_numbers[state_index]};
    }
    energies_t energy_of_state;
    if (has_energy()) {
        energy_of_state.emplace({energy_of_states->at(state_index)});
    }
    return derived().get_modified_copy(std::move(coefficients_of_state), std::move(energy_of_state),
                                       std::move(quantum_numbers_of_state));
}

template <typename Derived>
std::shared_ptr<const typename Basis<Derived>::ket_t>
Basis<Derived>::get_ket(size_t ket_index) const {
    return kets[ket_index];
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::begin() const {
    return kets.begin();
}

template <typename Derived>
typename Basis<Derived>::Iterator Basis<Derived>::end() const {
    return kets.end();
}

template <typename Derived>
Basis<Derived>::Iterator::Iterator(typename ketvec_t::const_iterator it) : it{std::move(it)} {}

template <typename Derived>
bool Basis<Derived>::Iterator::operator!=(const Iterator &other) const {
    return other.it != it;
}

template <typename Derived>
std::shared_ptr<const typename Basis<Derived>::ket_t> Basis<Derived>::Iterator::operator*() const {
    return *it;
}

template <typename Derived>
typename Basis<Derived>::Iterator &Basis<Derived>::Iterator::operator++() {
    ++it;
    return *this;
}

template <typename Derived>
size_t Basis<Derived>::get_number_of_states() const {
    return coefficients.cols();
}

template <typename Derived>
size_t Basis<Derived>::get_number_of_kets() const {
    return coefficients.rows();
}

template <typename Derived>
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
Basis<Derived>::get_sorter(const std::vector<SorterType> &labels) const {
    // Note that sorting by the energy is rejected by get_sorter_without_checks

    // Initialize the sorter
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> sorter(coefficients.cols());
    sorter.setIdentity();

    get_sorter_without_checks(labels, sorter);

    return sorter;
}

template <typename Derived>
std::vector<IndicesOfBlock>
Basis<Derived>::get_indices_of_blocks(const std::vector<SorterType> &labels) const {
    std::set<SorterType> unique_labels(labels.begin(), labels.end());
    perform_blocks_checks(unique_labels);

    // Get the blocks
    IndicesOfBlocksCreator blocks_creator({0, static_cast<size_t>(coefficients.cols())});
    get_indices_of_blocks_without_checks(unique_labels, blocks_creator);

    return blocks_creator.create();
}

template <typename Derived>
void Basis<Derived>::get_sorter_without_checks(
    const std::vector<SorterType> &labels,
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter) const {
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    // Currently, since states in a basis dont store the energy, they cannot be sorted by energy.
    // Checking this upfront also guarantees that only labels that refer to a quantum number are
    // encountered below.
    if (std::find(labels.begin(), labels.end(), SorterType::ENERGY) != labels.end()) {
        throw std::invalid_argument(
            "States in a basis do not store the energy and thus can not be sorted by it. "
            "Note that sorting a system by the energy is supported nevertheless.");
    }

    if (coefficients.cols() == 0) {
        return;
    }

    std::vector<const std::vector<real_t> *> quantum_numbers_of_labels;
    quantum_numbers_of_labels.reserve(labels.size());
    for (const auto &label : labels) {
        quantum_numbers_of_labels.push_back(&get_quantum_numbers_of_states(label));
    }

    int *perm_begin = sorter.indices().data();
    int *perm_end = perm_begin + coefficients.cols();

    // Sort the vector based on the requested labels
    set_task_status("Sorting basis states...");
    std::stable_sort(perm_begin, perm_end, [&](int a, int b) {
        for (const auto *quantum_numbers : quantum_numbers_of_labels) {
            if (std::abs((*quantum_numbers)[a] - (*quantum_numbers)[b]) > numerical_precision) {
                return (*quantum_numbers)[a] < (*quantum_numbers)[b];
            }
        }
        return false; // Elements are equal
    });
}

template <typename Derived>
void Basis<Derived>::get_indices_of_blocks_without_checks(
    const std::set<SorterType> &unique_labels, IndicesOfBlocksCreator &blocks_creator) const {
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    if (coefficients.cols() == 0) {
        return;
    }

    std::vector<const std::vector<real_t> *> quantum_numbers_of_labels;
    quantum_numbers_of_labels.reserve(unique_labels.size());
    for (const auto &label : unique_labels) {
        quantum_numbers_of_labels.push_back(&get_quantum_numbers_of_states(label));
    }

    // A new block starts whenever a state differs from its predecessor in one of the labels
    set_task_status("Identifying basis blocks...");
    for (int i = 1; i < coefficients.cols(); ++i) {
        for (const auto *quantum_numbers : quantum_numbers_of_labels) {
            if (std::abs((*quantum_numbers)[i] - (*quantum_numbers)[i - 1]) > numerical_precision) {
                blocks_creator.add(i);
                break;
            }
        }
    }
}

template <typename Derived>
std::shared_ptr<const Derived> Basis<Derived>::get_canonicalized_copy() const {
    auto [coefficients, energy_of_states, quantum_numbers_of_states] = get_canonical_states(kets);
    return derived().get_modified_copy(std::move(coefficients), std::move(energy_of_states),
                                       std::move(quantum_numbers_of_states));
}

template <typename Derived>
bool Basis<Derived>::is_canonical() const {
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    // The basis is canonical if and only if its coefficient matrix is the identity matrix, i.e.,
    // if the i-th state is the i-th ket
    if (coefficients.rows() != coefficients.cols()) {
        return false;
    }

    Eigen::Index num_ones = 0;
    for (Eigen::Index row = 0; row < coefficients.outerSize(); ++row) {
        for (typename Eigen::SparseMatrix<scalar_t, Eigen::RowMajor>::InnerIterator it(coefficients,
                                                                                       row);
             it; ++it) {
            if (it.row() == it.col() &&
                std::abs(it.value() - static_cast<scalar_t>(1)) <= numerical_precision) {
                ++num_ones;
            } else if (std::abs(it.value()) > numerical_precision) {
                return false;
            }
        }
    }
    return num_ones == coefficients.rows();
}

template <typename Derived>
std::shared_ptr<const Derived> Basis<Derived>::get_sorted_copy(
    const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter) const {
    // Create a copy of the current object
    auto result = std::make_shared<Derived>(derived());

    if (coefficients.cols() == 0) {
        return result;
    }

    // Apply the sorting
    set_task_status("Applying basis sorting...");
    result->coefficients = coefficients * sorter;

    set_task_status("Relabeling sorted basis states...");
    for (auto &[name, sorted_quantum_numbers] : result->quantum_numbers_of_states) {
        const auto &quantum_numbers = quantum_numbers_of_states.at(name);
        sorted_quantum_numbers.resize(sorter.size());
        for (int i = 0; i < sorter.size(); ++i) {
            sorted_quantum_numbers[i] = quantum_numbers[sorter.indices()[i]];
        }
    }
    if (has_energy()) {
        auto &sorted_energies = *result->energy_of_states;
        for (int i = 0; i < sorter.size(); ++i) {
            sorted_energies[i] = (*energy_of_states)[sorter.indices()[i]];
        }
    }

    return result;
}

// Explicit instantiations
template class Basis<BasisAtom<double>>;
template class Basis<BasisAtom<std::complex<double>>>;
template class Basis<BasisPair<double>>;
template class Basis<BasisPair<std::complex<double>>>;
} // namespace pairinteraction
