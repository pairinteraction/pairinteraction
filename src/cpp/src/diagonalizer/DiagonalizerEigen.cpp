#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>

namespace pairinteraction {

template <typename Scalar>
DiagonalizerEigen<Scalar>::DiagonalizerEigen(FPP fpp) : DiagonalizerInterface<Scalar>(fpp) {}

template <typename Scalar>
template <typename ScalarLim>
EigenSystemH<Scalar>
DiagonalizerEigen<Scalar>::dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                         double atol) const {
    using real_t = typename traits::NumTraits<Scalar>::real_t;

    // Subtract the mean of the diagonal elements from the diagonal
    real_t shift{};
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        this->template subtract_mean<ScalarLim>(matrix, shift, atol);

    // Diagonalize the shifted matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<ScalarLim>> eigensolver;
    eigensolver.compute(shifted_matrix);

    return {eigensolver.eigenvectors().sparseView(1, atol).template cast<Scalar>(),
            this->add_mean(eigensolver.eigenvalues(), shift)};
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerEigen<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                double atol) const {
    switch (this->fpp) {
    case FPP::FLOAT32:
        return dispatch_eigh<traits::restricted_t<Scalar, FPP::FLOAT32>>(matrix, atol);
    case FPP::FLOAT64:
        return dispatch_eigh<traits::restricted_t<Scalar, FPP::FLOAT64>>(matrix, atol);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

// Explicit instantiations
template class DiagonalizerEigen<double>;
template class DiagonalizerEigen<std::complex<double>>;
} // namespace pairinteraction
