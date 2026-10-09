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
#include <map>
#include <set>
#include <string>
#include <vector>

namespace pairinteraction {

template <typename Scalar>
class BasisAtom;

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
Basis<Derived>::Basis(ketvec_t &&kets)
    : kets(std::move(kets)), coefficients(static_cast<Eigen::Index>(this->kets.size()),
                                          static_cast<Eigen::Index>(this->kets.size())) {
    if (this->kets.empty()) {
        throw std::invalid_argument("The basis must contain at least one element.");
    }
    make_canonical();
}

template <typename Derived>
void Basis<Derived>::make_canonical() {
    auto n = static_cast<Eigen::Index>(kets.size());
    if (coefficients.rows() != n || coefficients.cols() != n) {
        throw std::invalid_argument("The coefficient matrix must be of shape n x n, where n is "
                                    "the number of kets.");
    }

    // Make the i-th state equal to the i-th ket
    coefficients.setIdentity();

    // Assign each state the quantum numbers of the ket it is equal to
    for (const auto &[label, name] : Derived::sorter_type_to_quantum_number_name) {
        std::vector<real_t> quantum_numbers;
        quantum_numbers.reserve(kets.size());
        for (const auto &ket : kets) {
            quantum_numbers.push_back(ket->has_quantum_number(name)
                                          ? static_cast<real_t>(ket->get_quantum_number(name))
                                          : std::numeric_limits<real_t>::max());
        }
        quantum_numbers_of_states[name] = std::move(quantum_numbers);
    }

    // Assign each state the energy of the ket it is equal to
    energy_of_states.clear();
    energy_of_states.reserve(kets.size());
    for (const auto &ket : kets) {
        energy_of_states.push_back(static_cast<real_t>(ket->get_energy()));
    }
}

template <typename Derived>
bool Basis<Derived>::has_quantum_number(const std::string &name) const {
    auto it = quantum_numbers_of_states.find(name);
    if (it == quantum_numbers_of_states.end()) {
        return false;
    }
    return std::none_of(it->second.begin(), it->second.end(), [](real_t quantum_number) {
        return quantum_number == std::numeric_limits<real_t>::max();
    });
}

template <typename Derived>
bool Basis<Derived>::has_energy() const {
    return std::none_of(energy_of_states.begin(), energy_of_states.end(),
                        [](real_t energy) { return energy == std::numeric_limits<real_t>::max(); });
}

template <typename Derived>
bool Basis<Derived>::is_sorted_by_energy() const {
    return has_energy() && std::is_sorted(energy_of_states.begin(), energy_of_states.end());
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
std::shared_ptr<const Derived> Basis<Derived>::copy_with_coefficients(
    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &values) const {
    if (values.rows() != coefficients.rows()) {
        throw std::invalid_argument("Incompatible number of rows.");
    }
    if (values.cols() != coefficients.cols()) {
        throw std::invalid_argument("Incompatible number of columns.");
    }

    // Create a copy of the current object and update the coefficients of the copy
    auto result = std::make_shared<Derived>(derived());

    result->coefficients = values;

    for (auto &[name, quantum_numbers] : result->quantum_numbers_of_states) {
        std::fill(quantum_numbers.begin(), quantum_numbers.end(),
                  std::numeric_limits<real_t>::max());
    }
    std::fill(result->energy_of_states.begin(), result->energy_of_states.end(),
              std::numeric_limits<real_t>::max());

    return result;
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_quantum_number(const std::string &name,
                                                                   size_t state_index) const {
    auto it = quantum_numbers_of_states.find(name);
    if (it == quantum_numbers_of_states.end()) {
        throw std::invalid_argument("The states are not labeled by the quantum number " + name +
                                    ".");
    }
    real_t quantum_number = it->second.at(state_index);
    if (quantum_number == std::numeric_limits<real_t>::max()) {
        throw std::invalid_argument("The state does not have a well-defined quantum number " +
                                    name + ".");
    }
    return quantum_number;
}

template <typename Derived>
typename Basis<Derived>::real_t Basis<Derived>::get_energy(size_t state_index) const {
    real_t energy = energy_of_states.at(state_index);
    if (energy == std::numeric_limits<real_t>::max()) {
        throw std::invalid_argument("The state does not have a well-defined energy.");
    }
    return energy;
}

template <typename Derived>
std::shared_ptr<const Derived> Basis<Derived>::get_state(size_t state_index) const {
    // Create a copy of the current object
    auto restricted = std::make_shared<Derived>(derived());

    // Restrict the copy to the single requested state
    restricted->coefficients = restricted->coefficients.col(state_index);

    for (auto &[name, quantum_numbers] : restricted->quantum_numbers_of_states) {
        quantum_numbers = {quantum_numbers_of_states.at(name)[state_index]};
    }
    restricted->energy_of_states = {energy_of_states.at(state_index)};

    return restricted;
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
std::shared_ptr<const Derived> Basis<Derived>::canonicalized() const {
    auto result = std::make_shared<Derived>(derived());
    size_t n = kets.size();
    result->coefficients.resize(n, n);
    result->make_canonical();
    return result;
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
std::shared_ptr<const Derived> Basis<Derived>::transformed(
    const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &sorter) const {
    // Create a copy of the current object
    auto transformed = std::make_shared<Derived>(derived());

    if (coefficients.cols() == 0) {
        return transformed;
    }

    // Apply the sorting
    set_task_status("Applying basis sorting...");
    transformed->coefficients = coefficients * sorter;

    set_task_status("Relabeling sorted basis states...");
    for (auto &[name, transformed_quantum_numbers] : transformed->quantum_numbers_of_states) {
        const auto &quantum_numbers = quantum_numbers_of_states.at(name);
        transformed_quantum_numbers.resize(sorter.size());
        for (int i = 0; i < sorter.size(); ++i) {
            transformed_quantum_numbers[i] = quantum_numbers[sorter.indices()[i]];
        }
    }
    transformed->energy_of_states.resize(sorter.size());
    for (int i = 0; i < sorter.size(); ++i) {
        transformed->energy_of_states[i] = energy_of_states[sorter.indices()[i]];
    }

    return transformed;
}

template <typename Derived>
std::shared_ptr<const Derived>
Basis<Derived>::transformed(const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &transformation,
                            const std::vector<real_t> &energy_of_transformed_states) const {
    if (!energy_of_transformed_states.empty() &&
        energy_of_transformed_states.size() != static_cast<size_t>(transformation.cols())) {
        throw std::invalid_argument(
            "The number of energies must match the number of transformed states.");
    }

    const real_t numerical_precision =
        100 * std::sqrt(coefficients.rows()) * std::numeric_limits<real_t>::epsilon();

    // Create a copy of the current object
    auto transformed = std::make_shared<Derived>(derived());

    if (coefficients.cols() == 0) {
        return transformed;
    }

    // Apply the transformation
    set_task_status("Applying basis transformation...");
    transformed->coefficients = coefficients * transformation;

    Eigen::SparseMatrix<real_t> probs = transformation.cwiseAbs2().transpose();
    // Drop explicitly stored zeros
    probs.prune(
        [](Eigen::Index /*row*/, Eigen::Index /*col*/, const real_t &value) { return value != 0; });
    Eigen::VectorX<real_t> norm = probs * Eigen::VectorX<real_t>::Ones(probs.cols());

    set_task_status("Updating transformed quantum numbers...");
    for (auto &[name, transformed_quantum_numbers] : transformed->quantum_numbers_of_states) {
        const auto &quantum_numbers = quantum_numbers_of_states.at(name);
        auto map = Eigen::Map<const Eigen::VectorX<real_t>>(
            quantum_numbers.data(), static_cast<Eigen::Index>(quantum_numbers.size()));
        Eigen::VectorX<real_t> val = probs * map;
        Eigen::VectorX<real_t> sq = probs * map.cwiseAbs2();
        transformed_quantum_numbers.resize(probs.rows());

        for (size_t i = 0; i < transformed_quantum_numbers.size(); ++i) {
            const real_t mean = val[i] / norm[i];
            const real_t mean_sq = sq[i] / norm[i];
            const real_t variance = std::abs(mean_sq - mean * mean);
            if (variance < numerical_precision * std::max(static_cast<real_t>(1), mean_sq)) {
                transformed_quantum_numbers[i] = std::round(mean * 2) / 2;
            } else {
                transformed_quantum_numbers[i] = std::numeric_limits<real_t>::max();
            }
        }
    }

    // The transformed states only have well-defined energies if they are explicitly specified,
    // e.g., as eigenenergies
    if (energy_of_transformed_states.empty()) {
        transformed->energy_of_states.assign(static_cast<size_t>(transformation.cols()),
                                             std::numeric_limits<real_t>::max());
    } else {
        transformed->energy_of_states = energy_of_transformed_states;
    }

    return transformed;
}

// Explicit instantiations
template class Basis<BasisAtom<double>>;
template class Basis<BasisAtom<std::complex<double>>>;
template class Basis<BasisPair<double>>;
template class Basis<BasisPair<std::complex<double>>>;
} // namespace pairinteraction
