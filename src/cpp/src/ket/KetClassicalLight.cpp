#include "pairinteraction/ket/KetClassicalLight.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <fmt/core.h>
#include <memory>
#include <string>

namespace pairinteraction {
template <typename Real>
KetClassicalLight<Real>::KetClassicalLight(Private /*unused*/, Real photon_energy, int q)
    : Ket<Real>(photon_energy * q, 0, 0, Parity::ODD), photon_energy(photon_energy),
      quantum_number_q(q) {}

template <typename Real>
Real KetClassicalLight<Real>::get_photon_energy() const {
    return photon_energy;
}

template <typename Real>
std::string KetClassicalLight<Real>::get_label() const {
    std::string label;
    label += fmt::format("{:d}", quantum_number_q);
    label += fmt::format(",{:g}GHz", photon_energy);
    return label;
}

template <typename Real>
std::shared_ptr<KetClassicalLight<Real>>
KetClassicalLight<Real>::get_ket_for_different_quantum_number_m(Real new_quantum_number_m) const {
    if (new_quantum_number_m != 0) {
        throw std::invalid_argument(
            "Classical light cannot have a different quantum number m than zero.");
    }
    return std::make_shared<KetClassicalLight<Real>>(*this);
}

template <typename Real>
int KetClassicalLight<Real>::get_quantum_number_q() const {
    return quantum_number_q;
}

template <typename Real>
bool KetClassicalLight<Real>::operator==(const KetClassicalLight<Real> &other) const {
    return Ket<Real>::operator==(other) && photon_energy == other.photon_energy &&
        quantum_number_q == other.quantum_number_q;
}

template <typename Real>
bool KetClassicalLight<Real>::operator!=(const KetClassicalLight<Real> &other) const {
    return !(*this == other);
}

template <typename Real>
size_t KetClassicalLight<Real>::hash::operator()(const KetClassicalLight<Real> &k) const {
    size_t seed = typename Ket<Real>::hash()(k);
    utils::hash_combine(seed, k.photon_energy);
    utils::hash_combine(seed, k.quantum_number_q);
    return seed;
}

// Explicit instantiations
template class KetClassicalLight<double>;
} // namespace pairinteraction
