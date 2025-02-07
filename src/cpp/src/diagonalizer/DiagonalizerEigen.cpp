#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace pairinteraction {

template <typename Scalar, typename ScalarRstr>
EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                   int precision) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    const real_t shift = matrix.diagonal().real().mean();
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> identity(dim, dim);
    identity.setIdentity();
    Eigen::MatrixX<ScalarRstr> shifted_matrix =
        (matrix - shift * identity).template cast<ScalarRstr>();

    // Diagonalize the shifted matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<ScalarRstr>> eigensolver;
    eigensolver.compute(shifted_matrix);

    // Add the mean of the diagonal elements back to the eigenvalues
    Eigen::VectorX<real_t> eigenvalues = eigensolver.eigenvalues().template cast<real_t>();
    eigenvalues.array() += shift;

    return {
        eigensolver.eigenvectors().sparseView(std::pow(10, -precision), 1).template cast<Scalar>(),
        eigenvalues};
}

template <typename Scalar>
DiagonalizerEigen<Scalar>::DiagonalizerEigen(FPP fpp) : DiagonalizerInterface<Scalar>(fpp) {}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerEigen<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                int precision) const {
    switch (this->fpp) {
    case FPP::FLOAT32:
        return dispatch_eigh<Scalar, traits::restricted_t<Scalar, FPP::FLOAT32>>(matrix, precision);
    case FPP::FLOAT64:
        return dispatch_eigh<Scalar, traits::restricted_t<Scalar, FPP::FLOAT64>>(matrix, precision);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

// Explicit instantiations
template class DiagonalizerEigen<double>;
template class DiagonalizerEigen<std::complex<double>>;
} // namespace pairinteraction
