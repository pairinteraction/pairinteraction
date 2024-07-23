#include "ket/KetClassicalLightCreator.hpp"

#include "ket/KetClassicalLight.hpp"

template <typename Real>
KetClassicalLightCreator<Real>::KetClassicalLightCreator(Real photon_energy, int q)
    : photon_energy(photon_energy), quantum_number_q(q) {}

template <typename Real>
KetClassicalLightCreator<Real> &KetClassicalLightCreator<Real>::set_photon_energy(Real value) {
    photon_energy.emplace(value);
    return *this;
}

template <typename Real>
KetClassicalLightCreator<Real> &KetClassicalLightCreator<Real>::set_quantum_number_q(int value) {
    quantum_number_q.emplace(value);
    return *this;
}

template <typename Real>
std::shared_ptr<const KetClassicalLight<Real>> KetClassicalLightCreator<Real>::create() const {
    return std::make_shared<KetClassicalLight<Real>>(
        typename KetClassicalLight<Real>::Private(),
        photon_energy.value_or(std::numeric_limits<Real>::quiet_NaN()),
        quantum_number_q.value_or(std::numeric_limits<int>::max()));
}

// Explicit instantiations
template class KetClassicalLightCreator<float>;
template class KetClassicalLightCreator<double>;
