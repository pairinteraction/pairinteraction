#pragma once

#include "interfaces/TransformationBuilderInterface.hpp"
#include "utils/eigen_assertion.hpp"
#include "utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <array>
#include <complex>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

enum class TransformationType : unsigned char;

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
class Basis
    : public TransformationBuilderInterface<typename traits::CrtpTraits<Derived>::scalar_t> {
public:
    using scalar_t = typename traits::CrtpTraits<Derived>::scalar_t;
    using real_t = typename traits::CrtpTraits<Derived>::real_t;
    using ket_t = typename traits::CrtpTraits<Derived>::ket_t;
    using ketvec_t = typename traits::CrtpTraits<Derived>::ketvec_t;

    Basis() = delete;
    virtual ~Basis() = default;

    const ketvec_t &get_kets() const;

    size_t get_number_of_states() const;
    size_t get_number_of_kets() const;
    real_t get_quantum_number_f(size_t index_state) const;
    real_t get_quantum_number_m(size_t index_state) const;
    int get_parity(size_t index_state) const;
    std::shared_ptr<const ket_t> get_ket_with_largest_overlap(size_t index_state) const;
    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_coefficients() const;

    class Iterator {
    public:
        Iterator(typename ketvec_t::const_iterator it);
        bool operator!=(const Iterator &other) const;
        std::shared_ptr<const ket_t> operator*() const;
        Iterator &operator++();

    private:
        typename ketvec_t::const_iterator it;
    };

    Iterator begin() const;
    Iterator end() const;

    const Transformation<scalar_t> &get_transformation() const override;
    Transformation<scalar_t> get_rotator(real_t alpha, real_t beta, real_t gamma) const override;
    Sorting get_sorter(TransformationType label) const override;
    Blocks get_blocks(TransformationType label) const override;

    Sorting get_sorter_without_checks(TransformationType label) const;
    Blocks get_blocks_without_checks(TransformationType label) const;

    std::shared_ptr<Derived> transformed(const Transformation<scalar_t> &transformation) const;
    std::shared_ptr<Derived> transformed(const Sorting &transformation) const;

protected:
    Basis(ketvec_t &&kets);
    ketvec_t kets;
    Transformation<scalar_t> coefficients;
    std::unordered_map<size_t, size_t> ket_id_to_index;

private:
    const Derived &derived() const;
    std::vector<real_t> quantum_number_f_of_states;
    std::vector<real_t> quantum_number_m_of_states;
    std::vector<int> parity_of_states;
    std::vector<int> ket_of_states;
};
