#include "pairinteraction/diagonalizer/DiagonalizerEigen.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace pairinteraction {
template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerEigen<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                int precision) const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<Scalar>> eigensolver;
    eigensolver.compute(matrix);

    return {eigensolver.eigenvectors().sparseView(std::pow(10, -precision), 1),
            eigensolver.eigenvalues()};
}

// Explicit instantiations
template class DiagonalizerEigen<double>;
template class DiagonalizerEigen<std::complex<double>>;
} // namespace pairinteraction
