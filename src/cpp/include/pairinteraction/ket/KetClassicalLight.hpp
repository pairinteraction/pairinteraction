#pragma once

#include "pairinteraction/ket/Ket.hpp"

#include <string>
#include <type_traits>

namespace pairinteraction {
template <typename Real>
class KetClassicalLightCreator;

template <typename Real>
class KetClassicalLight : public Ket<Real> {
    static_assert(std::is_floating_point_v<Real>);

    friend class KetClassicalLightCreator<Real>;
    struct Private {};

public:
    KetClassicalLight(Private /*unused*/, Real photon_energy, int q);

    std::string get_label() const override;
    std::shared_ptr<KetClassicalLight<Real>>
    get_ket_for_different_quantum_number_m(Real new_quantum_number_m) const;
    Real get_photon_energy() const;
    int get_quantum_number_q() const;

    bool operator==(const KetClassicalLight<Real> &other) const;

    struct hash {
        std::size_t operator()(const KetClassicalLight<Real> &k) const;
    };

private:
    Real photon_energy;
    int quantum_number_q;
};

extern template class KetClassicalLight<float>;
extern template class KetClassicalLight<double>;
} // namespace pairinteraction
