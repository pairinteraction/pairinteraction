#pragma once

#include <optional>

#include "ket/KetClassicalLight.hpp"

template <typename Real>
class KetClassicalLight;

template <typename Real>
class KetClassicalLightCreator {
public:
    KetClassicalLightCreator();
    KetClassicalLightCreator(Real energy, int n_sigma_p, int n_pi, int n_sigma_m);
    KetClassicalLightCreator<Real> &set_energy(Real value);
    KetClassicalLightCreator<Real> &set_quantum_number_n_sigma_p(int value);
    KetClassicalLightCreator<Real> &set_quantum_number_n_pi(int value);
    KetClassicalLightCreator<Real> &set_quantum_number_n_sigma_m(int value);
    KetClassicalLight<Real> create() const;

private:
    std::optional<Real> energy;
    std::optional<int> quantum_number_n_sigma_p;
    std::optional<int> quantum_number_n_pi;
    std::optional<int> quantum_number_n_sigma_m;
    std::optional<int> parity;
};

extern template class KetClassicalLightCreator<float>;
extern template class KetClassicalLightCreator<double>;
