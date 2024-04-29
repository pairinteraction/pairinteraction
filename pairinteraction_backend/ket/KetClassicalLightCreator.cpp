#include "ket/KetClassicalLightCreator.hpp"

#include "ket/KetClassicalLight.hpp"

template <typename Real>
KetClassicalLightCreator<Real>::KetClassicalLightCreator(Real energy, int q)
    : energy(energy), quantum_number_q(q) {}

template <typename Real>
KetClassicalLightCreator<Real> &KetClassicalLightCreator<Real>::set_energy(Real value) {
    energy.emplace(value);
    return *this;
}

template <typename Real>
KetClassicalLightCreator<Real> &KetClassicalLightCreator<Real>::set_quantum_number_q(int value) {
    quantum_number_q.emplace(value);
    return *this;
}

template <typename Real>
KetClassicalLight<Real> KetClassicalLightCreator<Real>::create() const {
    return KetClassicalLight<Real>(energy.value_or(std::numeric_limits<Real>::quiet_NaN()),
                                   quantum_number_q.value_or(std::numeric_limits<int>::max()),
                                   1000); // TODO use a unique id instead 1000
}

template class KetClassicalLightCreator<float>;
template class KetClassicalLightCreator<double>;
