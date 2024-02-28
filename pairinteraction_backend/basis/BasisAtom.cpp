#include "BasisAtom.hpp"
#include "ket/KetAtom.hpp"

template <typename T, bool is_complex>
BasisAtom<T, is_complex>::BasisAtom() : Basis<T, is_complex>() {}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_energy(T min, T max) {
    this->ensure_not_assembled();
    min_energy = min;
    max_energy = max;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_quantum_number_f(float min, float max) {
    this->ensure_not_assembled();
    min_quantum_number_f = min;
    max_quantum_number_f = max;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_quantum_number_m(float min, float max) {
    this->ensure_not_assembled();
    min_quantum_number_m = min;
    max_quantum_number_m = max;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_parity(int parity) {
    this->ensure_not_assembled();
    this->parity = parity;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_quantum_number_n(int n) {
    this->ensure_not_assembled();
    quantum_number_n = n;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_quantum_number_nu(T min, T max) {
    this->ensure_not_assembled();
    min_quantum_number_nu = min;
    max_quantum_number_nu = max;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_quantum_number_l(T min, T max) {
    this->ensure_not_assembled();
    min_quantum_number_l = min;
    max_quantum_number_l = max;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_quantum_number_s(T min, T max) {
    this->ensure_not_assembled();
    min_quantum_number_s = min;
    max_quantum_number_s = max;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::restrict_quantum_number_j(T min, T max) {
    this->ensure_not_assembled();
    min_quantum_number_j = min;
    max_quantum_number_j = max;
}

template <typename T, bool is_complex>
void BasisAtom<T, is_complex>::ensure_assembled_kets() {
    // TODO perform database request
}

// Explicit instantiations
template class BasisAtom<float, false>;
template class BasisAtom<double, false>;
template class BasisAtom<float, true>;
template class BasisAtom<double, true>;
