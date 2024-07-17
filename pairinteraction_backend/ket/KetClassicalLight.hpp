#pragma once

#include <string>
#include <type_traits>

#include "ket/Ket.hpp"

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
    size_t get_id() const override;
    size_t get_id_for_different_quantum_number_m(Real new_quantum_number_m) const override;
    Real get_photon_energy() const;
    int get_quantum_number_q() const;

private:
    Real photon_energy;
    int quantum_number_q;
};

extern template class KetClassicalLight<float>;
extern template class KetClassicalLight<double>;
