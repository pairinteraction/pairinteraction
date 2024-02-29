#pragma once

#include <string>

#include "ket/Ket.hpp"
#include "ket/KetAtomCreator.hpp"

template <typename T>
class KetAtomCreator;

/**
 * @class KetAtom
 *
 * @brief Class for representing atomic kets.
 *
 * @tparam T Real number type.
 */
template <typename T>
class KetAtom : public Ket<T> {
public:
    std::string get_species() const;
    int get_quantum_number_n() const;
    T get_quantum_number_nu() const;
    T get_quantum_number_l() const;
    T get_quantum_number_s() const;
    T get_quantum_number_j() const;

private:
    friend class KetAtomCreator<T>;
    KetAtom(T energy, float f, float m, int p, std::string label, std::string species, int n,
            T nu_exp, T nu_std, T l_exp, T l_std, T s_exp, T s_std, T j_exp, T j_std);
    std::string species;
    int quantum_number_n;
    T quantum_number_nu_exp;
    T quantum_number_nu_std;
    T quantum_number_l_exp;
    T quantum_number_l_std;
    T quantum_number_s_exp;
    T quantum_number_s_std;
    T quantum_number_j_exp;
    T quantum_number_j_std;
};

extern template class KetAtom<float>;
extern template class KetAtom<double>;
