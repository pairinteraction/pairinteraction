#ifndef KETATOM_HPP
#define KETATOM_HPP

#include <limits>
#include <string>

#include "Ket.hpp"

template <typename T>
class KetAtom : public Ket<T> {
public:
    KetAtom(std::string species);
    KetAtom(T energy, float f, float m, int p, std::string label, std::string species, int n,
            T nu_exp, T nu_std, T l_exp, T l_std, T s_exp, T s_std, T j_exp, T j_std);
    int get_quantum_number_n();
    T get_quantum_number_nu();
    T get_quantum_number_l();
    T get_quantum_number_s();
    T get_quantum_number_j();
    void set_energy(T value);
    void set_quantum_number_f(float value);
    void set_quantum_number_m(float value);
    void set_parity(int value);
    void set_quantum_number_n(int value);
    void set_quantum_number_nu(T value);
    void set_quantum_number_l(T value);
    void set_quantum_number_s(T value);
    void set_quantum_number_j(T value);

private:
    void ensure_assembled_ket() override;
    std::string species{""};
    int quantum_number_n{std::numeric_limits<int>::max()};
    T quantum_number_nu_exp{std::numeric_limits<T>::max()};
    T quantum_number_nu_std{std::numeric_limits<T>::max()};
    T quantum_number_l_exp{std::numeric_limits<T>::max()};
    T quantum_number_l_std{std::numeric_limits<T>::max()};
    T quantum_number_s_exp{std::numeric_limits<T>::max()};
    T quantum_number_s_std{std::numeric_limits<T>::max()};
    T quantum_number_j_exp{std::numeric_limits<T>::max()};
    T quantum_number_j_std{std::numeric_limits<T>::max()};
};

extern template class KetAtom<float>;
extern template class KetAtom<double>;

#endif // KETATOM_HPP
