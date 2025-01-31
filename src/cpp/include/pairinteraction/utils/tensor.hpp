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
calculate_tensor_product(const std::shared_ptr<const BasisPair<Scalar>> &basis,
                         const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix1,
                         const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix2);

extern template Eigen::SparseMatrix<float, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<float>> &,
                         const Eigen::SparseMatrix<float, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<float, Eigen::RowMajor> &);
extern template Eigen::SparseMatrix<double, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<double>> &,
                         const Eigen::SparseMatrix<double, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<double, Eigen::RowMajor> &);
extern template Eigen::SparseMatrix<std::complex<float>, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<std::complex<float>>> &,
                         const Eigen::SparseMatrix<std::complex<float>, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<std::complex<float>, Eigen::RowMajor> &);
extern template Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor>
calculate_tensor_product(const std::shared_ptr<const BasisPair<std::complex<double>>> &,
                         const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> &,
                         const Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> &);
} // namespace pairinteraction::utils
