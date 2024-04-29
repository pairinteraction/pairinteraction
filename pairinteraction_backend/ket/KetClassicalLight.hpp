#pragma once

#include <string>

#include "ket/Ket.hpp"

template <typename Real>
class KetClassicalLightCreator;

template <typename Real>
class KetClassicalLight : public Ket<Real> {
public:
    std::string get_label() const override;
    size_t get_id() const override;
    size_t get_id_for_different_quantum_number_m(float new_quantum_number_m) const override;
    Real get_photon_energy() const;
    int get_quantum_number_q() const;

private:
    friend class KetClassicalLightCreator<Real>;
    KetClassicalLight(Real photon_energy, int q, size_t id);
    Real photon_energy;
    size_t id;
    int quantum_number_q;
};

extern template class KetClassicalLight<float>;
extern template class KetClassicalLight<double>;
