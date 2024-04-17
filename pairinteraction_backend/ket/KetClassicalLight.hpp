#pragma once

#include <string>

#include "KetClassicalLightCreator.hpp"
#include "ket/Ket.hpp"
#include "ket/KetClassicalLightCreator.hpp"

template <typename Real>
class KetClassicalLightCreator;

template <typename Real>
class KetClassicalLight : public Ket<Real> {
public:
    int get_quantum_number_n_sigma_p() const;
    int get_quantum_number_n_pi() const;
    int get_quantum_number_n_sigma_m() const;

private:
    friend class KetClassicalLightCreator<Real>;
    KetClassicalLight(Real energy, int n_sigma_p, int n_pi, int n_sigma_m, std::string label,
                      size_t id); // p can only be -1 for photons
                                  // But this I handle in KetClassicalLightCreator.
    int quantum_number_n_sigma_p;
    int quantum_number_n_pi;
    int quantum_number_n_sigma_m;
};

extern template class KetClassicalLight<float>;
extern template class KetClassicalLight<double>;
