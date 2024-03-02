#pragma once

#include <Eigen/SparseCore>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "ket/Ket.hpp"
#include "utils/Traits.hpp"

/**
 * @class Basis
 *
 * @brief Base class for a basis
 *
 * This base class represents a basis. It comprises a list of ket states and a matrix of
 * coefficients. Using CRPT, it is a base class for specific basis implementations. Its
 * constructor is protected to indicate that derived classes should not allow direct instantiation.
 * Instead, a factory class should be provided that is a friend of the derived class and can create
 * instances of it.
 *
 * @tparam Derived Derived class.
 */

template <typename Derived>
class Basis {
public:
    using scalar_t = typename internal::BasisTraits<Derived>::scalar_t;
    using real_t = typename internal::BasisTraits<Derived>::real_t;
    using ket_t = typename internal::BasisTraits<Derived>::ket_t;
    using ketvec_t = typename internal::BasisTraits<Derived>::ketvec_t;

    Basis() = delete;
    const ket_t &get_ket(size_t index) const;
    real_t get_energy(size_t index) const;
    float get_quantum_number_f(size_t index) const;
    float get_quantum_number_m(size_t index) const;
    int get_parity(size_t index) const;
    std::string get_label(size_t index) const;

    class Iterator {
    public:
        Iterator(const Basis<Derived> &basis, size_t index);
        bool operator!=(const Iterator &other) const;
        const ket_t &operator*() const;
        Iterator &operator++();

    private:
        const Basis<Derived> &basis;
        size_t index;
    };

    Iterator begin() const;
    Iterator end() const;

protected:
    Basis(ketvec_t &&kets);

private:
    const Derived &derived() const;
    std::vector<real_t> energies;
    std::vector<float> quantum_numbers_f;
    std::vector<float> quantum_numbers_m;
    std::vector<int> parities;
    std::vector<std::string> labels;
    Eigen::SparseMatrix<scalar_t> coefficients;
    bool is_standard_basis;
    ketvec_t kets;
};
