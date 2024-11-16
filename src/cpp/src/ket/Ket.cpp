#include "pairinteraction/ket/Ket.hpp"

#include "pairinteraction/enums/Parity.hpp"
#include "pairinteraction/utils/hash.hpp"

#include <limits>

namespace pairinteraction {
template <typename Real>
Ket<Real>::Ket(Real energy, Real f, Real m, Parity p)
    : energy(energy), quantum_number_f(f), quantum_number_m(m), parity(p) {}

template <typename Real>
bool Ket<Real>::has_quantum_number_f() const {
    return quantum_number_f != std::numeric_limits<Real>::max();
}

template <typename Real>
bool Ket<Real>::has_quantum_number_m() const {
    return quantum_number_m != std::numeric_limits<Real>::max();
}

template <typename Real>
bool Ket<Real>::has_parity() const {
    return parity != Parity::UNKNOWN;
}

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

template <typename Real>
bool Ket<Real>::operator==(const Ket<Real> &other) const {
    return energy == other.energy && quantum_number_f == other.quantum_number_f &&
        quantum_number_m == other.quantum_number_m && parity == other.parity;
}

template <typename Real>
size_t Ket<Real>::hash::operator()(const Ket<Real> &k) const {
    size_t seed = 0;
    utils::hash_combine(seed, k.energy);
    utils::hash_combine(seed, k.quantum_number_f);
    utils::hash_combine(seed, k.quantum_number_m);
    utils::hash_combine(seed, static_cast<int>(k.parity));
    return seed;
}

// Explicit instantiations
template class Ket<float>;
template class Ket<double>;
} // namespace pairinteraction
