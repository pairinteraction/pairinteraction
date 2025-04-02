// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/ket/KetPair.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <limits>
#include <string>

namespace pairinteraction {
template <typename Scalar>
KetPair<Scalar>::KetPair(
    Private /*unused*/, std::initializer_list<size_t> atomic_indices,
    std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases, real_t energy)
    : Ket(energy, calculate_quantum_number_f(atomic_indices, atomic_bases),
          calculate_quantum_number_m(atomic_indices, atomic_bases),
          calculate_parity(atomic_indices, atomic_bases)),
      atomic_indices(atomic_indices), atomic_bases(atomic_bases) {
    if (atomic_indices.size() != atomic_bases.size()) {
        throw std::invalid_argument(
            "The number of atomic indices, and atomic bases must be the same.");
    }
}

template <typename Scalar>
std::string KetPair<Scalar>::get_label() const {
    std::string label;
    std::string separator;
    for (size_t atom_index = 0; atom_index < atomic_indices.size(); ++atom_index) {
        label += separator +
            atomic_bases[atom_index]
                ->get_corresponding_ket(atomic_indices[atom_index])
                ->get_label();
        separator = "; ";
    }
    return label;
}

template <typename Scalar>
std::shared_ptr<KetPair<Scalar>>
KetPair<Scalar>::get_ket_for_different_quantum_number_m(real_t /*new_quantum_number_m*/) const {
    // If we use symmetrized states so that the quantum_number_f is the total
    // angular quantum number, the quantum_number_m is the magnetic quantum number
    // corresponding to the total angular quantum number and we can implement this
    // method.
    throw std::runtime_error("Not implemented.");
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
    return Ket::operator==(other) && atomic_indices == other.atomic_indices &&
        atomic_bases == other.atomic_bases;
}

template <typename Scalar>
bool KetPair<Scalar>::operator!=(const KetPair<Scalar> &other) const {
    return !(*this == other);
}

template <typename Scalar>
size_t KetPair<Scalar>::hash::operator()(const KetPair<Scalar> &k) const {
    size_t seed = typename Ket::hash()(k);
    for (const auto &index : k.atomic_indices) {
        utils::hash_combine(seed, index);
    }
    for (const auto &basis : k.atomic_bases) {
        utils::hash_combine(seed, reinterpret_cast<std::uintptr_t>(basis.get()));
    }
    return seed;
}

template <typename Scalar>
typename KetPair<Scalar>::real_t KetPair<Scalar>::calculate_quantum_number_f(
    const std::vector<size_t> & /*indices*/,
    const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> & /*bases*/) {
    // Because this ket state is not symmetrized, the quantum_number_f is not well-defined.
    return std::numeric_limits<real_t>::max();
}

template <typename Scalar>
typename KetPair<Scalar>::real_t KetPair<Scalar>::calculate_quantum_number_m(
    const std::vector<size_t> &indices,
    const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> &bases) {
    for (const auto &basis : bases) {
        if (!basis->has_quantum_number_m()) {
            return std::numeric_limits<real_t>::max();
        }
    }
    real_t total_quantum_number_m = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        total_quantum_number_m += bases[i]->get_quantum_number_m(indices[i]);
    }
    return total_quantum_number_m;
}

template <typename Scalar>
Parity KetPair<Scalar>::calculate_parity(
    const std::vector<size_t> & /*indices*/,
    const std::vector<std::shared_ptr<const BasisAtom<Scalar>>> & /*bases*/) {
    // Because this ket state is not symmetrized, the parity is not well-defined.
    return Parity::UNKNOWN;
}

// Explicit instantiations
template class KetPair<double>;
template class KetPair<std::complex<double>>;
} // namespace pairinteraction
