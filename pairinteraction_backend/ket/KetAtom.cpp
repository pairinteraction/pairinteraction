#include "ket/KetAtom.hpp"
#include <string>

template <typename Real>
KetAtom<Real>::KetAtom(Real energy, float f, float m, int p, std::string label, size_t id,
                       std::string species, int n, Real nu_exp, Real nu_std, Real l_exp, Real l_std,
                       Real s_exp, Real s_std, Real j_exp, Real j_std)
    : Ket<Real>(energy, f, m, p, label, id), species(species), quantum_number_n(n),
      quantum_number_nu_exp(nu_exp), quantum_number_nu_std(nu_std), quantum_number_l_exp(l_exp),
      quantum_number_l_std(l_std), quantum_number_s_exp(s_exp), quantum_number_s_std(s_std),
      quantum_number_j_exp(j_exp), quantum_number_j_std(j_std) {}

template <typename Real>
std::string KetAtom<Real>::get_species() const {
    return species;
}

template <typename Real>
int KetAtom<Real>::get_quantum_number_n() const {
    return quantum_number_n;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_nu() const {
    return quantum_number_nu_exp;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_l() const {
    return quantum_number_l_exp;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_s() const {
    return quantum_number_s_exp;
}

template <typename Real>
Real KetAtom<Real>::get_quantum_number_j() const {
    return quantum_number_j_exp;
}

// Explicit instantiations
template class KetAtom<float>;
template class KetAtom<double>;
