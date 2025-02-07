#include "pairinteraction/ket/KetClassicalLight.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <fmt/core.h>
#include <memory>
#include <string>

namespace pairinteraction {
KetClassicalLight::KetClassicalLight(Private /*unused*/, double photon_energy, int q)
    : Ket(photon_energy * q, 0, 0, Parity::ODD), photon_energy(photon_energy), quantum_number_q(q) {
}

double KetClassicalLight::get_photon_energy() const { return photon_energy; }

std::string KetClassicalLight::get_label() const {
    std::string label;
    label += fmt::format("{:d}", quantum_number_q);
    label += fmt::format(",{:g}GHz", photon_energy);
    return label;
}

std::shared_ptr<KetClassicalLight>
KetClassicalLight::get_ket_for_different_quantum_number_m(double new_quantum_number_m) const {
    if (new_quantum_number_m != 0) {
        throw std::invalid_argument(
            "Classical light cannot have a different quantum number m than zero.");
    }
    return std::make_shared<KetClassicalLight>(*this);
}

int KetClassicalLight::get_quantum_number_q() const { return quantum_number_q; }

bool KetClassicalLight::operator==(const KetClassicalLight &other) const {
    return Ket::operator==(other) && photon_energy == other.photon_energy &&
        quantum_number_q == other.quantum_number_q;
}

bool KetClassicalLight::operator!=(const KetClassicalLight &other) const {
    return !(*this == other);
}

size_t KetClassicalLight::hash::operator()(const KetClassicalLight &k) const {
    size_t seed = typename Ket::hash()(k);
    utils::hash_combine(seed, k.photon_energy);
    utils::hash_combine(seed, k.quantum_number_q);
    return seed;
}
} // namespace pairinteraction
