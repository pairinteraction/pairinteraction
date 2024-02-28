#pragma once

#include <Eigen/SparseCore>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "ket/Ket.hpp"

template <typename T>
struct make_real {
    using type = T;
};

template <typename T>
struct make_real<std::complex<T>> {
    using type = T;
};

template <typename T>
using real_t = typename make_real<T>::type;

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
    Basis(std::vector<std::shared_ptr<const Ket<real_t<T>>>> &&kets);
    std::vector<std::shared_ptr<const Ket<real_t<T>>>> kets;

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
