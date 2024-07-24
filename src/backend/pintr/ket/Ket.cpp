#include "pintr/ket/Ket.hpp"

#include "pintr/enums/Parity.hpp"

#include <limits>

template <typename Real>
Ket<Real>::Ket(Real energy, Real f, Real m, Parity p)
    : energy(energy), quantum_number_f(f), quantum_number_m(m), parity(p) {}

template <typename Real>
Real Ket<Real>::get_energy() const {
    return energy;
}

template <typename Real>
Real Ket<Real>::get_quantum_number_f() const {
    return quantum_number_f;
}

template <typename Real>
Real Ket<Real>::get_quantum_number_m() const {
    return quantum_number_m;
}

template <typename Real>
Parity Ket<Real>::get_parity() const {
    return parity;
}

// Explicit instantiations
template class Ket<float>;
template class Ket<double>;
