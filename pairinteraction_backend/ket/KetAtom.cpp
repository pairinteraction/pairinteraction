#include "ket/KetAtom.hpp"
#include <string>

template <typename T>
KetAtom<T>::KetAtom(T energy, float f, float m, int p, std::string label, std::string species,
                    int n, T nu_exp, T nu_std, T l_exp, T l_std, T s_exp, T s_std, T j_exp, T j_std)
    : Ket<T>(energy, f, m, p, label), species(species), quantum_number_n(n),
      quantum_number_nu_exp(nu_exp), quantum_number_nu_std(nu_std), quantum_number_l_exp(l_exp),
      quantum_number_l_std(l_std), quantum_number_s_exp(s_exp), quantum_number_s_std(s_std),
      quantum_number_j_exp(j_exp), quantum_number_j_std(j_std) {}

template <typename T>
std::string KetAtom<T>::get_species() const {
    return species;
}

template <typename T>
int KetAtom<T>::get_quantum_number_n() const {
    return quantum_number_n;
}

template <typename T>
T KetAtom<T>::get_quantum_number_nu() const {
    return quantum_number_nu_exp;
}

template <typename T>
T KetAtom<T>::get_quantum_number_l() const {
    return quantum_number_l_exp;
}

template <typename T>
T KetAtom<T>::get_quantum_number_s() const {
    return quantum_number_s_exp;
}

template <typename T>
T KetAtom<T>::get_quantum_number_j() const {
    return quantum_number_j_exp;
}

// Explicit instantiations
template class KetAtom<float>;
template class KetAtom<double>;
