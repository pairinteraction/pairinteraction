#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <array>
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
    using scalar_t = typename Traits::BasisTraits<Derived>::scalar_t;
    using real_t = typename Traits::BasisTraits<Derived>::real_t;
    using ket_t = typename Traits::BasisTraits<Derived>::ket_t;
    using ketvec_t = typename Traits::BasisTraits<Derived>::ketvec_t;

    Basis() = delete;
    const ket_t &get_ket(size_t index_ket) const;
    real_t get_energy(size_t index_state) const;

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
    std::vector<int> get_sorter_according_to_kets() const;
    std::vector<int> get_sorter_according_to_energies() const;
    std::vector<int> get_sorter_according_to_quantum_number_m() const;
    std::vector<int> get_sorter_according_to_parity() const;

    void transform(const Eigen::SparseMatrix<scalar_t> &transformator);
    void rotate(real_t alpha, real_t beta, real_t gamma);
    void rotate(std::array<real_t, 3> to_z_axis, std::array<real_t, 3> to_y_axis);
    void sort(const std::vector<int> &sorter);
    void sort_according_to_kets();
    void sort_according_to_energies();
    void sort_according_to_quantum_number_m();
    void sort_according_to_parity();

    void set_eigen_basis(Eigen::SparseMatrix<scalar_t> evecs, std::vector<real_t> evals);

protected:
    Basis(ketvec_t &&kets);

private:
    const Derived &derived() const;
    std::vector<real_t> energy_of_states;
    Eigen::SparseMatrix<scalar_t> coefficients;
    bool is_standard_basis;
    ketvec_t kets;
};
