#include "ket/KetClassicalLight.hpp"

#include <fmt/format.h>
#include <limits>
#include <string>

#include "KetClassicalLight.hpp"
#include "utils/hash.hpp"

template <typename Real>
KetClassicalLight<Real>::KetClassicalLight(Real photon_energy, int q)
    : Ket<Real>(photon_energy * q, 0, 0, -1), photon_energy(photon_energy), quantum_number_q(q) {}

template <typename Real>
Real KetClassicalLight<Real>::get_photon_energy() const {
    return photon_energy;
}

template <typename Real>
std::string KetClassicalLight<Real>::get_label() const {
    std::string label = "";
    label += fmt::format("{:d}", quantum_number_q);
    label += fmt::format("_{{{:g} GHz}}", photon_energy);
    return label;
}

template <typename Real>
size_t KetClassicalLight<Real>::get_id() const {
    size_t seed = 0;
    hash::hash_combine(seed, quantum_number_q);
    hash::hash_combine(seed, photon_energy);
    return seed;
}

template <typename Real>
size_t
KetClassicalLight<Real>::get_id_for_different_quantum_number_m(float new_quantum_number_m) const {
    if (new_quantum_number_m != 0) {
        throw std::invalid_argument(
            "Classical light cannot have a different quantum number m than zero.");
    }
    return this->get_id();
}

template <typename Real>
int KetClassicalLight<Real>::get_quantum_number_q() const {
    return quantum_number_q;
}

template class KetClassicalLight<float>;
template class KetClassicalLight<double>;
