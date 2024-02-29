#pragma once

#include "basis/BasisAtom.hpp"
#include <complex>
#include <optional>
#include <string>

template <typename T>
class BasisAtom;

/**
 * @class BasisAtomCreator
 *
 * @brief Builder class for creating BasisAtom objects.
 *
 * @tparam T Complex number type.
 */
template <typename T>
class BasisAtomCreator {
public:
    BasisAtomCreator(std::string species);
    void restrict_energy(T min, T max);
    void restrict_quantum_number_f(float min, float max);
    void restrict_quantum_number_m(float min, float max);
    void restrict_parity(int parity);
    void restrict_quantum_number_n(int min, int max);
    void restrict_quantum_number_nu(T min, T max);
    void restrict_quantum_number_l(T min, T max);
    void restrict_quantum_number_s(T min, T max);
    void restrict_quantum_number_j(T min, T max);
    BasisAtom<T> create() const;

private:
    std::string species;
    std::optional<T> min_energy;
    std::optional<T> max_energy;
    std::optional<float> min_quantum_number_f;
    std::optional<float> max_quantum_number_f;
    std::optional<float> min_quantum_number_m;
    std::optional<float> max_quantum_number_m;
    std::optional<int> parity;
    std::optional<int> min_quantum_number_n;
    std::optional<int> max_quantum_number_n;
    std::optional<T> min_quantum_number_nu;
    std::optional<T> max_quantum_number_nu;
    std::optional<T> min_quantum_number_l;
    std::optional<T> max_quantum_number_l;
    std::optional<T> min_quantum_number_s;
    std::optional<T> max_quantum_number_s;
    std::optional<T> min_quantum_number_j;
    std::optional<T> max_quantum_number_j;
};

extern template class BasisAtomCreator<float>;
extern template class BasisAtomCreator<double>;
extern template class BasisAtomCreator<std::complex<float>>;
extern template class BasisAtomCreator<std::complex<double>>;
