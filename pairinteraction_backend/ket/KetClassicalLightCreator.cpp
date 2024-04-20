#include "ket/KetClassicalLightCreator.hpp"

#include <exception>
#include <limits>
#include <string>

template <typename Real>
KetClassicalLightCreator<Real>::KetClassicalLightCreator() : parity(-1) {}

template <typename Real>
KetClassicalLightCreator<Real>::KetClassicalLightCreator(Real energy, int n_sigma_p, int n_pi,
                                                         int n_sigma_m)
    : energy(energy), quantum_number_n_sigma_p(n_sigma_p), quantum_number_n_pi(n_pi),
      quantum_number_n_sigma_m(n_sigma_m), parity(-1) {}

template <typename Real>
KetClassicalLightCreator<Real> &KetClassicalLightCreator<Real>::set_energy(Real value) {
    energy.emplace(value);
    return *this;
}

template <typename Real>
KetClassicalLightCreator<Real> &
KetClassicalLightCreator<Real>::set_quantum_number_n_sigma_p(int value) {
    quantum_number_n_sigma_p.emplace(value);
    return *this;
}

template <typename Real>
KetClassicalLightCreator<Real> &KetClassicalLightCreator<Real>::set_quantum_number_n_pi(int value) {
    quantum_number_n_pi.emplace(value);
    return *this;
}

template <typename Real>
KetClassicalLightCreator<Real> &
KetClassicalLightCreator<Real>::set_quantum_number_n_sigma_m(int value) {
    quantum_number_n_sigma_m.emplace(value);
    return *this;
}

template <typename Real>
KetClassicalLight<Real> KetClassicalLightCreator<Real>::create() const {
    return KetClassicalLight<Real>(
        energy.value_or(std::numeric_limits<Real>::quiet_NaN()),
        quantum_number_n_sigma_p.value_or(std::numeric_limits<int>::max()),
        quantum_number_n_pi.value_or(std::numeric_limits<int>::max()),
        quantum_number_n_sigma_m.value_or(std::numeric_limits<int>::max()),
        1000); // TODO use a unique id instead 1000
}

template class KetClassicalLightCreator<float>;
template class KetClassicalLightCreator<double>;
