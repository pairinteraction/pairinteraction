#pragma once

#include <limits>
#include <optional>
#include <string>

#include "ket/KetAtom.hpp"

template <typename T>
class KetAtom;

/**
 * @class KetAtomCreator
 *
 * @brief Builder class for creating KetAtom objects.
 *
 * @tparam T Real number type.
 */
template <typename T>
class KetAtomCreator {
public:
    KetAtomCreator(std::string species);
    KetAtomCreator(std::string species, int n, T l, float j, float m);
    void set_energy(T value);
    void set_quantum_number_f(float value);
    void set_quantum_number_m(float value);
    void set_parity(int value);
    void set_quantum_number_n(int value);
    void set_quantum_number_nu(T value);
    void set_quantum_number_l(T value);
    void set_quantum_number_s(T value);
    void set_quantum_number_j(T value);
    KetAtom<T> create() const;

private:
    std::string species;
    std::optional<T> energy;
    std::optional<float> quantum_number_f;
    std::optional<float> quantum_number_m;
    std::optional<int> parity;
    std::optional<int> quantum_number_n;
    std::optional<T> quantum_number_nu;
    std::optional<T> quantum_number_l;
    std::optional<T> quantum_number_s;
    std::optional<T> quantum_number_j;
};

extern template class KetAtomCreator<float>;
extern template class KetAtomCreator<double>;
