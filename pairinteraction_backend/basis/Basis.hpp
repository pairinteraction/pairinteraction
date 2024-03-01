#pragma once

#include <Eigen/SparseCore>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "ket/Ket.hpp"
#include "utils/Traits.hpp"

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

// TODO: CRTP

// https://pixorblog.wordpress.com/2019/08/14/curiously-recurring-template-pattern-crtp-in-depth/

// https://www.fluentcpp.com/2020/09/11/replacing-crtp-static-polymorphism-with-concepts/

// https://pixorblog.wordpress.com/2019/08/14/curiously-recurring-template-pattern-crtp-in-depth/

// https://stackoverflow.com/questions/76229523/why-is-a-nested-c-template-argument-not-deduced-even-though-it-is-constrained

#include "ket/KetAtom.hpp"

template <typename Derived> class Basis;

namespace internal {

template <typename Derived>
struct traits<Basis<Derived>> : traits<Derived> {};

}

// template <typename Derived>
// class Basis;

// template <typename Scalar>
// template <class T> requires std::derived_from<T, Basis<T<Scalar>>>
// struct DerivedTraits
// {
//     using Scalar = Scalar;
//     using Real = real_t<Scalar>;
//     using Ket = Ket<Real>;
// };

template <typename Derived>
class Basis {
public:
    using Scalar = typename internal::traits<Derived>::Scalar;   // TODO remove
    using Real = real_t<Scalar>;                                 // TODO remove
    using KetType = typename internal::traits<Derived>::KetType; // TODO remove
    //using KetType = decltype(std::declval<Derived>().get_ket(std::declval<size_t>())); // TODO remove

    const Derived &derived() const;
    const auto &get_ket(size_t index) const;
    Real get_energy(size_t index) const;
    float get_quantum_number_f(size_t index) const;
    float get_quantum_number_m(size_t index) const;
    int get_parity(size_t index) const;
    std::string get_label(size_t index) const;

    class Iterator {
    public:
        Iterator(const Basis<Derived> &basis, size_t index);
        bool operator!=(const Iterator &other) const;
        const KetType &operator*() const; // TODO auto
        Iterator &operator++();

    private:
        const Basis<Derived> &basis;
        size_t index;
    };

    Iterator begin() const;
    Iterator end() const;

protected:
    Basis();

private:
    std::vector<Real> energies;
    std::vector<float> quantum_numbers_f;
    std::vector<float> quantum_numbers_m;
    std::vector<int> parities;
    std::vector<std::string> labels;
    Eigen::SparseMatrix<Scalar> coefficients;
    bool is_standard_basis;
};
