// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <string>

namespace pairinteraction {
template <typename Scalar>
KetPair<Scalar>::KetPair(
    Private /*unused*/, std::initializer_list<size_t> atomic_indices,
    std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases, real_t energy,
    std::unordered_map<std::string, double> quantum_numbers)
    : Ket(energy), quantum_numbers(std::move(quantum_numbers)), atomic_indices(atomic_indices),
      atomic_bases(atomic_bases) {
    if (atomic_indices.size() != atomic_bases.size()) {
        throw std::invalid_argument(
            "The number of atomic indices, and atomic bases must be the same.");
    }
}

template <typename Scalar>
bool KetPair<Scalar>::has_quantum_number(const std::string &name) const {
    return quantum_numbers.contains(name);
}

template <typename Scalar>
double KetPair<Scalar>::get_quantum_number(const std::string &name) const {
    return quantum_numbers.at(name);
}

template <typename Scalar>
std::vector<std::shared_ptr<const BasisAtom<Scalar>>> KetPair<Scalar>::get_atomic_states() const {
    std::vector<std::shared_ptr<const BasisAtom<Scalar>>> atomic_states;
    atomic_states.reserve(atomic_indices.size());
    for (size_t atom_index = 0; atom_index < atomic_indices.size(); ++atom_index) {
        atomic_states.push_back(atomic_bases[atom_index]->get_state(atomic_indices[atom_index]));
    }
    return atomic_states;
}

template <typename Scalar>
bool KetPair<Scalar>::operator==(const KetPair<Scalar> &other) const {
    return Ket::operator==(other) && quantum_numbers == other.quantum_numbers &&
        atomic_indices == other.atomic_indices && atomic_bases == other.atomic_bases;
}

template <typename Scalar>
bool KetPair<Scalar>::operator!=(const KetPair<Scalar> &other) const {
    return !(*this == other);
}

template <typename Scalar>
size_t KetPair<Scalar>::hash::operator()(const KetPair<Scalar> &k) const {
    size_t seed = typename Ket::hash()(k);
    // The quantum numbers are stored in an unordered map, so we combine the per-entry hashes in an
    // order-independent way (via xor) to obtain a deterministic result.
    size_t quantum_numbers_hash = 0;
    for (const auto &[key, value] : k.quantum_numbers) {
        size_t entry_seed = 0;
        utils::hash_combine(entry_seed, key);
        utils::hash_combine(entry_seed, value);
        quantum_numbers_hash ^= entry_seed;
    }
    utils::hash_combine(seed, quantum_numbers_hash);
    for (const auto &index : k.atomic_indices) {
        utils::hash_combine(seed, index);
    }
    for (const auto &basis : k.atomic_bases) {
        utils::hash_combine(seed, reinterpret_cast<std::uintptr_t>(basis.get()));
    }
    return seed;
}

// Explicit instantiations
template class KetPair<double>;
template class KetPair<std::complex<double>>;
} // namespace pairinteraction
