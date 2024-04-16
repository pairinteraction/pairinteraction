#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <array>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include "ket/Ket.hpp"
#include "utils/traits.hpp"

/**
 * @class Basis
 *
 * @brief Base class for a basis
 *
 * This base class represents a basis. It comprises a list of ket states and a matrix of
 * coefficients. The rows of the coefficient matrix correspond to indices of ket states and
 * the columns to indices of basis vectors.
 * Using CRPT, it is a base class for specific basis implementations. Its
 * constructor is protected to indicate that derived classes should not allow direct instantiation.
 * Instead, a factory class should be provided that is a friend of the derived class and can create
 * instances of it.
 *
 * @tparam Derived Derived class.
 */

template <typename Derived>
class Basis {
public:
    using scalar_t = typename traits::BasisTraits<Derived>::scalar_t;
    using real_t = typename traits::BasisTraits<Derived>::real_t;
    using ket_t = typename traits::BasisTraits<Derived>::ket_t;
    using ketvec_t = typename traits::BasisTraits<Derived>::ketvec_t;

    Basis() = delete;
    size_t get_number_of_states() const;
    size_t get_number_of_kets() const;
    const ket_t &get_ket(size_t index_ket) const;
    float get_quantum_number_f(size_t index_state) const;
    float get_quantum_number_m(size_t index_state) const;
    int get_parity(size_t index_state) const;

    enum class Label : unsigned char {
        NONE = 1 << 0,
        QUANTUM_NUMBER_F = 1 << 1,
        QUANTUM_NUMBER_M = 1 << 2,
        PARITY = 1 << 3,
        KET = 1 << 4
    };

    friend inline constexpr Label operator&(Label x, Label y) {
        return static_cast<Label>(static_cast<unsigned char>(x) & static_cast<unsigned char>(y));
    }

    friend inline constexpr Label operator|(Label x, Label y) {
        return static_cast<Label>(static_cast<unsigned char>(x) | static_cast<unsigned char>(y));
    }

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

    Eigen::SparseMatrix<scalar_t> get_rotator(real_t alpha, real_t beta, real_t gamma) const;
    Eigen::SparseMatrix<scalar_t> get_rotator(std::array<real_t, 3> to_z_axis,
                                              std::array<real_t, 3> to_y_axis) const;
    std::vector<int> get_sorter(Label label) const;
    std::vector<int> get_blocks(Label label) const;

    void transform(const Eigen::SparseMatrix<scalar_t> &transformator);
    void rotate(real_t alpha, real_t beta, real_t gamma);
    void rotate(std::array<real_t, 3> to_z_axis, std::array<real_t, 3> to_y_axis);
    void sort(const std::vector<int> &sorter);
    void sort(Label label);

protected:
    Basis(ketvec_t &&kets);

private:
    const Derived &derived() const;
    std::vector<float> quantum_number_f_of_states;
    std::vector<float> quantum_number_m_of_states;
    std::vector<int> parity_of_states;
    std::vector<int> ket_of_states;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> coefficients;
    ketvec_t kets;
    bool is_standard_basis;
    Label sortation;
};
