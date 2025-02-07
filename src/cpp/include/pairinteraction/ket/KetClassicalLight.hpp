#pragma once

#include "pairinteraction/ket/Ket.hpp"

#include <string>
#include <type_traits>

namespace pairinteraction {
class KetClassicalLightCreator;

class KetClassicalLight : public Ket {
    friend class KetClassicalLightCreator;
    struct Private {};

public:
    KetClassicalLight(Private /*unused*/, double photon_energy, int q);

    std::string get_label() const override;
    std::shared_ptr<KetClassicalLight>
    get_ket_for_different_quantum_number_m(double new_quantum_number_m) const;
    double get_photon_energy() const;
    int get_quantum_number_q() const;

    bool operator==(const KetClassicalLight &other) const;
    bool operator!=(const KetClassicalLight &other) const;

    struct hash {
        std::size_t operator()(const KetClassicalLight &k) const;
    };

private:
    double photon_energy;
    int quantum_number_q;
};

} // namespace pairinteraction
