#pragma once

#include <Eigen/SparseCore>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "ket/Ket.hpp"

/**
 * @struct make_real
 *
 * @brief Helper struct to extract the real type from a scalar type.
 *
 * @tparam T The type from which to extract the real type.
 */

template <typename T>
struct make_real {
    using type = T;
};

template <typename T>
struct make_real<std::complex<T>> {
    using type = T;
};

/**
 * @typedef real_t
 *
 * @brief Helper type to extract the real type from a scalar type
 *
 * @tparam T The type from which to extract the real type.
 */

template <typename T>
using real_t = typename make_real<T>::type;

/**
 * @class Basis
 *
 * @brief Base class for a basis
 *
 * This base class represents a basis. It comprises a list of ket states and a matrix of
 * coefficients. It is meant to be used as a base class for specific basis implementations. Its
 * constructor is protected to indicate that derived classes should not allow direct instantiation.
 * Instead, a factory class should be provided that is a friend of the derived class and can create
 * instances of it.
 *
 * @tparam T Complex number type.
 */

template <typename T>
class Basis {
public:
    Basis() = delete;
    const Ket<real_t<T>> &get_ket(size_t index) const;
    real_t<T> get_energy(size_t index) const;
    float get_quantum_number_f(size_t index) const;
    float get_quantum_number_m(size_t index) const;
    int get_parity(size_t index) const;
    std::string get_label(size_t index) const;

    class Iterator {
    public:
        Iterator(const Basis<T> &basis, size_t index);
        bool operator!=(const Iterator &other) const;
        const Ket<real_t<T>> &operator*() const;
        Iterator &operator++();

    private:
        const Basis<T> &basis;
        size_t index;
    };

    Iterator begin() const;
    Iterator end() const;

protected:
    using KetPtrVec = std::vector<std::shared_ptr<const Ket<real_t<T>>>>;
    Basis(KetPtrVec &&kets);
    KetPtrVec kets;

private:
    std::vector<real_t<T>> energies;
    std::vector<float> quantum_numbers_f;
    std::vector<float> quantum_numbers_m;
    std::vector<int> parities;
    std::vector<std::string> labels;
    Eigen::SparseMatrix<T> coefficients;
    bool is_standard_basis;
};

extern template class Basis<float>;
extern template class Basis<double>;
extern template class Basis<std::complex<float>>;
extern template class Basis<std::complex<double>>;
