#include "ket/KetClassicalLight.hpp"
#include <limits>

template <typename Real>
KetClassicalLight<Real>::KetClassicalLight(Real energy, int n_sigma_p, int n_pi, int n_sigma_m,
                                           size_t id)
    : Ket<Real>(energy, std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), -1,
                id),
      quantum_number_n_sigma_p(n_sigma_p), quantum_number_n_pi(n_pi),
      quantum_number_n_sigma_m(n_sigma_m) {}

template <typename Real>
std::string KetClassicalLight<Real>::get_label() const {
    return "TODO";
}

template <typename Real>
int KetClassicalLight<Real>::get_quantum_number_n_sigma_p() const {
    return quantum_number_n_sigma_p;
}

template <typename Real>
int KetClassicalLight<Real>::get_quantum_number_n_pi() const {
    return quantum_number_n_pi;
}

template <typename Real>
int KetClassicalLight<Real>::get_quantum_number_n_sigma_m() const {
    return quantum_number_n_sigma_m;
}

template class KetClassicalLight<float>;
template class KetClassicalLight<double>;
