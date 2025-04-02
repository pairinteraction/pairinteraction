// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <memory>

namespace pairinteraction {
template <typename Scalar>
class BasisPair;
} // namespace pairinteraction

namespace pairinteraction::utils {
template <typename Scalar>
Eigen::SparseMatrix<Scalar, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<Scalar>> &basis_initial,
                         const std::shared_ptr<const BasisPair<Scalar>> &basis_final,
                         const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix1,
                         const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix2);

extern template Eigen::SparseMatrix<double, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<double>> &,
                         const std::shared_ptr<const BasisPair<double>> &,
                         const Eigen::SparseMatrix<double, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<double, Eigen::RowMajor> &);
extern template Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<std::complex<double>>> &,
                         const std::shared_ptr<const BasisPair<std::complex<double>>> &,
                         const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> &);
} // namespace pairinteraction::utils
