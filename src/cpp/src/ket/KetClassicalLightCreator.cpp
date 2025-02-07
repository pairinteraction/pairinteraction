#include "pairinteraction/ket/KetClassicalLightCreator.hpp"

#include "pairinteraction/ket/KetClassicalLight.hpp"

#include <limits>

namespace pairinteraction {
KetClassicalLightCreator::KetClassicalLightCreator(double photon_energy, int q)
    : photon_energy(photon_energy), quantum_number_q(q) {}

KetClassicalLightCreator &KetClassicalLightCreator::set_photon_energy(double value) {
    photon_energy.emplace(value);
    return *this;
}

KetClassicalLightCreator &KetClassicalLightCreator::set_quantum_number_q(int value) {
    quantum_number_q.emplace(value);
    return *this;
}

std::shared_ptr<const KetClassicalLight> KetClassicalLightCreator::create() const {
    return std::make_shared<KetClassicalLight>(
        typename KetClassicalLight::Private(),
        photon_energy.value_or(std::numeric_limits<double>::quiet_NaN()),
        quantum_number_q.value_or(std::numeric_limits<int>::max()));
}
} // namespace pairinteraction
