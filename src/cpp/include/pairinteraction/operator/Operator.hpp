// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <memory>

namespace pairinteraction {
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
Derived &operator+=(Operator<Derived> &lhs, const Operator<Derived> &rhs);

template <typename Derived>
Derived &operator-=(Operator<Derived> &lhs, const Operator<Derived> &rhs);

template <typename Derived>
Derived operator+(const Operator<Derived> &lhs, const Operator<Derived> &rhs);

template <typename Derived>
Derived operator-(const Operator<Derived> &lhs, const Operator<Derived> &rhs);

template <typename Derived>
class Operator
    : public TransformationBuilderInterface<typename traits::CrtpTraits<Derived>::scalar_t> {
public:
    using scalar_t = typename traits::CrtpTraits<Derived>::scalar_t;
    using real_t = typename traits::CrtpTraits<Derived>::real_t;
    using ketvec_t = typename traits::CrtpTraits<Derived>::ketvec_t;
    using basis_t = typename traits::CrtpTraits<Derived>::basis_t;

    Operator(std::shared_ptr<const basis_t> basis);
    virtual ~Operator() = default;

    std::shared_ptr<const basis_t> get_basis() const;
    std::shared_ptr<const basis_t> &get_basis();

    const Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_matrix() const;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &get_matrix();

    const Transformation<scalar_t> &get_transformation() const override;
    Transformation<scalar_t> get_rotator(real_t alpha, real_t beta, real_t gamma) const override;
    Sorting get_sorter(const std::vector<TransformationType> &labels) const override;
    std::vector<IndicesOfBlock>
    get_indices_of_blocks(const std::vector<TransformationType> &labels) const override;

    void get_indices_of_blocks_without_checks(const std::set<TransformationType> &unique_labels,
                                              IndicesOfBlocksCreator &blocks) const;

    Derived transformed(const Transformation<scalar_t> &transformation) const;
    Derived transformed(const Sorting &transformation) const;

    friend Derived operator*
        <>(const typename Operator<Derived>::scalar_t &lhs, const Operator<Derived> &rhs);
    friend Derived operator*
        <>(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs);
    friend Derived operator/
        <>(const Operator<Derived> &lhs, const typename Operator<Derived>::scalar_t &rhs);
    // clang-format off
    friend Derived &operator+=<>(Operator<Derived> &lhs, const Operator<Derived> &rhs);
    friend Derived &operator-=<>(Operator<Derived> &lhs, const Operator<Derived> &rhs);
    friend Derived operator+<>(const Operator<Derived> &lhs, const Operator<Derived> &rhs);
    friend Derived operator-<>(const Operator<Derived> &lhs, const Operator<Derived> &rhs);
    // clang-format on

protected:
    void initialize_as_energy_operator();
    void initialize_from_matrix(Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> &&matrix);
    std::shared_ptr<const basis_t> basis;
    Eigen::SparseMatrix<scalar_t, Eigen::RowMajor> matrix;

private:
    const Derived &derived() const;
    Derived &derived_mutable();
};
} // namespace pairinteraction
