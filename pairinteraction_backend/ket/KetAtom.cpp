#include "KetAtom.hpp"

template <typename T>
KetAtom<T>::KetAtom(std::string species) : species(species) {}

template <typename T>
KetAtom<T>::KetAtom(T energy, float f, float m, int p, std::string label, std::string species,
                    int n, T nu_exp, T nu_std, T l_exp, T l_std, T s_exp, T s_std, T j_exp, T j_std)
    : Ket<T>(energy, f, m, p, label), species(species), quantum_number_n(n),
      quantum_number_nu_exp(nu_exp), quantum_number_nu_std(nu_std), quantum_number_l_exp(l_exp),
      quantum_number_l_std(l_std), quantum_number_s_exp(s_exp), quantum_number_s_std(s_std),
      quantum_number_j_exp(j_exp), quantum_number_j_std(j_std) {}

template <typename T>
void KetAtom<T>::ensure_assembled_ket() {
    // TODO perform database request
}

template <typename T>
int KetAtom<T>::get_quantum_number_n() {
    this->ensure_assembled();
    return quantum_number_n;
}

template <typename T>
T KetAtom<T>::get_quantum_number_nu() {
    this->ensure_assembled();
    return quantum_number_nu_exp;
}

template <typename T>
T KetAtom<T>::get_quantum_number_l() {
    this->ensure_assembled();
    return quantum_number_l_exp;
}

template <typename T>
T KetAtom<T>::get_quantum_number_s() {
    this->ensure_assembled();
    return quantum_number_s_exp;
}

template <typename T>
T KetAtom<T>::get_quantum_number_j() {
    this->ensure_assembled();
    return quantum_number_j_exp;
}

template <typename T>
void KetAtom<T>::set_energy(T value) {
    this->ensure_not_assembled();
    this->energy = value;
}

template <typename T>
void KetAtom<T>::set_quantum_number_f(float value) {
    this->ensure_not_assembled();
    this->quantum_number_f = value;
}

template <typename T>
void KetAtom<T>::set_quantum_number_m(float value) {
    this->ensure_not_assembled();
    this->quantum_number_m = value;
}

template <typename T>
void KetAtom<T>::set_parity(int value) {
    this->ensure_not_assembled();
    this->parity = value;
}

template <typename T>
void KetAtom<T>::set_quantum_number_n(int value) {
    this->ensure_not_assembled();
    this->quantum_number_n = value;
}

template <typename T>
void KetAtom<T>::set_quantum_number_nu(T value) {
    this->ensure_not_assembled();
    this->quantum_number_nu_exp = value;
}

template <typename T>
void KetAtom<T>::set_quantum_number_l(T value) {
    this->ensure_not_assembled();
    this->quantum_number_l_exp = value;
}

template <typename T>
void KetAtom<T>::set_quantum_number_s(T value) {
    this->ensure_not_assembled();
    this->quantum_number_s_exp = value;
}

template <typename T>
void KetAtom<T>::set_quantum_number_j(T value) {
    this->ensure_not_assembled();
    this->quantum_number_j_exp = value;
}

// Explicit instantiations
template class KetAtom<float>;
template class KetAtom<double>;
