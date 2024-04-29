#include "ket/KetClassicalLight.hpp"
#include "KetClassicalLight.hpp"
#include <limits>

template <typename Real>
KetClassicalLight<Real>::KetClassicalLight(Real photon_energy, int q, size_t id)
    : Ket<Real>(photon_energy * q, 0, 0, -1), id(id), photon_energy(photon_energy),
      quantum_number_q(q) {}

template <typename Real>
Real KetClassicalLight<Real>::get_photon_energy() const {
    return photon_energy;
}

template <typename Real>
std::string KetClassicalLight<Real>::get_label() const {
    return "TODO";
}

template <typename Real>
size_t KetClassicalLight<Real>::get_id() const {
    return id;
}

template <typename Real>
size_t
KetClassicalLight<Real>::get_id_for_different_quantum_number_m(float new_quantum_number_m) const {
    if (new_quantum_number_m != 0) {
        throw std::invalid_argument(
            "Classical light cannot have a different quantum number m than zero.");
    }
    return id;
}

template <typename Real>
int KetClassicalLight<Real>::get_quantum_number_q() const {
    return quantum_number_q;
}

template class KetClassicalLight<float>;
template class KetClassicalLight<double>;
