#include "pairinteraction/ket/KetCombined.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <limits>
#include <string>

namespace pairinteraction {
template <typename Scalar>
KetCombined<Scalar>::KetCombined(
    Private /*unused*/, std::initializer_list<size_t> atomic_indices,
    std::initializer_list<std::string> atomic_labels,
    std::initializer_list<std::shared_ptr<const BasisAtom<Scalar>>> atomic_bases, real_t energy,
    real_t quantum_number_f, real_t quantum_number_m, Parity parity)
    : Ket<real_t>(energy, quantum_number_f, quantum_number_m, parity),
      atomic_indices(atomic_indices), atomic_labels(atomic_labels), atomic_bases(atomic_bases) {
    if (atomic_indices.size() != atomic_labels.size() ||
        atomic_indices.size() != atomic_bases.size()) {
        throw std::invalid_argument(
            "The number of atomic indices, atomic labels, and atomic bases must be the same.");
    }
}

template <typename Scalar>
std::string KetCombined<Scalar>::get_label() const {
    std::string label;
    std::string separator;
    for (const auto &atomic_label : atomic_labels) {
        label += separator + atomic_label;
        separator = "; ";
    }
    return label;
}

template <typename Scalar>
size_t KetCombined<Scalar>::get_id() const {
    if (atomic_indices.empty()) {
        return 0;
    }
    size_t linear_index = atomic_indices[0];
    for (size_t k = 1; k < atomic_indices.size(); ++k) {
        linear_index = linear_index * atomic_bases[k]->get_number_of_states() + atomic_indices[k];
    }
    return linear_index;
}

template <typename Scalar>
size_t
KetCombined<Scalar>::get_id_for_different_quantum_number_m(real_t /*new_quantum_number_m*/) const {
    // If we use symmetrized states so that the quantum_number_f is the total
    // angular quantum number, the quantum_number_m is the magnetic quantum number
    // corresponding to the total angular quantum number and we can implement this
    // method.
    throw std::runtime_error("Not implemented.");
}

// Explicit instantiations
template class KetCombined<float>;
template class KetCombined<double>;
template class KetCombined<std::complex<float>>;
template class KetCombined<std::complex<double>>;
} // namespace pairinteraction
