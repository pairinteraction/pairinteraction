#include "pairinteraction/ket/KetCombined.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <limits>
#include <string>

namespace pairinteraction {
template <typename Real>
KetCombined<Real>::KetCombined(
    Private /*unused*/, size_t id, Real energy, Real quantum_number_f, Real quantum_number_m,
    Parity parity, std::vector<std::shared_ptr<const Ket<Real>>> &&kets_with_largest_overlap)
    : Ket<Real>(energy, quantum_number_f, quantum_number_m, parity), id(id),
      kets_with_largest_overlap(std::move(kets_with_largest_overlap)) {}

template <typename Real>
std::string KetCombined<Real>::get_label() const {
    std::string label;
    std::string separator;
    for (const auto &ket : kets_with_largest_overlap) {
        label += separator + ket->get_label();
        separator = "; ";
    }
    return label;
}

template <typename Real>
size_t KetCombined<Real>::get_id() const {
    return id;
}

template <typename Real>
size_t
KetCombined<Real>::get_id_for_different_quantum_number_m(Real /*new_quantum_number_m*/) const {
    // If we use symmetrized states so that the quantum_number_f is the total
    // angular quantum number, the quantum_number_m is the magnetic quantum number
    // corresponding to the total angular quantum number and we can implement this
    // method.
    throw std::runtime_error("Not implemented.");
}

// Explicit instantiations
template class KetCombined<float>;
template class KetCombined<double>;
} // namespace pairinteraction
