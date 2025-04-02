// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/diagonalize/DiagonalizerEigen.hpp"

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>

namespace pairinteraction {

template <typename Scalar>
DiagonalizerEigen<Scalar>::DiagonalizerEigen(FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type) {}

template <typename Scalar>
template <typename ScalarLim>
EigenSystemH<Scalar>
DiagonalizerEigen<Scalar>::dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                         double rtol) const {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    real_t shift{};
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        this->template subtract_mean<ScalarLim>(matrix, shift, rtol);

    // Diagonalize the shifted matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<ScalarLim>> eigensolver;
    eigensolver.compute(shifted_matrix);

    return {eigensolver.eigenvectors()
                .sparseView(1, 0.5 * rtol / std::sqrt(dim))
                .template cast<Scalar>(),
            this->add_mean(eigensolver.eigenvalues(), shift)};
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerEigen<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                double rtol) const {
    switch (this->float_type) {
    case FloatType::FLOAT32:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT32>>(matrix, rtol);
    case FloatType::FLOAT64:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT64>>(matrix, rtol);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

// Explicit instantiations
template class DiagonalizerEigen<double>;
template class DiagonalizerEigen<std::complex<double>>;
} // namespace pairinteraction
