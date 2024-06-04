#pragma once

#include <Eigen/SparseCore>
#include <array>
#include <memory>
#include <vector>

#include "interfaces/TransformationBuilderInterface.hpp"
#include "utils/traits.hpp"

enum class TransformationType : unsigned char;

template <typename Derived>
class Operator;

template <typename Derived>
Derived operator*(const typename Operator<Derived>::scalar_t &lhs, const Operator<Derived> &rhs);

template <typename Derived>
Derived operator*(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs);

template <typename Derived>
Derived operator/(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs);

template <typename Derived>
Derived operator+(const Operator<Derived> &lhs, const Operator<Derived> &rhs);

template <typename Derived>
Derived operator-(const Operator<Derived> &lhs, const Operator<Derived> &rhs);

template <typename Derived>
class Operator
    : public TransformationBuilderInterface<typename traits::OperatorTraits<Derived>::scalar_t> {
public:
    using scalar_t = typename traits::OperatorTraits<Derived>::scalar_t;
    using real_t = typename traits::OperatorTraits<Derived>::real_t;
    using ketvec_t = typename traits::OperatorTraits<Derived>::ketvec_t;
    using basis_t = typename traits::OperatorTraits<Derived>::basis_t;

    Operator(std::shared_ptr<const basis_t> basis);

    std::shared_ptr<const basis_t> get_basis() const;

    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_matrix() const;

    const Transformation<scalar_t> &get_transformation() const override;
    Transformation<scalar_t> get_rotator(real_t alpha, real_t beta, real_t gamma) const override;
    Sorting get_sorter(TransformationType label) const override;
    Blocks get_blocks(TransformationType label) const override;

    Derived transform(const Transformation<scalar_t> &transformation) const;
    Derived transform(const Sorting &transformation) const;

    friend Derived operator*
        <>(const typename Operator<Derived>::scalar_t &lhs, const Operator<Derived> &rhs);
    friend Derived operator*
        <>(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs);
    friend Derived operator/
        <>(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs);
    friend Derived operator+<>(const Operator<Derived> &lhs, const Operator<Derived> &rhs);
    friend Derived operator-<>(const Operator<Derived> &lhs, const Operator<Derived> &rhs);

protected:
    std::shared_ptr<const basis_t> basis;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> matrix;

private:
    const Derived &derived() const;
};
