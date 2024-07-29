#include "pintr/diagonalizer/DiagonalizerEigen.hpp"

#include "pintr/utils/eigen_assertion.hpp"
#include "pintr/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace pintr {
template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerEigen<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                real_t min_eigenvalue, real_t max_eigenvalue, int precision) const {
    int dim = matrix.rows();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<Scalar>> eigensolver;
    eigensolver.compute(matrix);

    Eigen::VectorX<real_t> evals = eigensolver.eigenvalues();
    Eigen::MatrixX<Scalar> evecs = eigensolver.eigenvectors();

    if (min_eigenvalue > std::numeric_limits<real_t>::lowest() / 2 ||
        max_eigenvalue < std::numeric_limits<real_t>::max() / 2) {
        auto *it_begin = std::lower_bound(evals.data(), evals.data() + dim, min_eigenvalue);
        auto *it_end = std::upper_bound(evals.data(), evals.data() + dim, max_eigenvalue);
        evecs = evecs
                    .block(0, std::distance(evals.data(), it_begin), dim,
                           std::distance(it_begin, it_end))
                    .eval();
        evals =
            evals.segment(std::distance(evals.data(), it_begin), std::distance(it_begin, it_end))
                .eval();
    }

    return {evecs.sparseView(std::pow(10, -precision), 1), evals};
}

// Explicit instantiations
template class DiagonalizerEigen<float>;
template class DiagonalizerEigen<double>;
template class DiagonalizerEigen<std::complex<float>>;
template class DiagonalizerEigen<std::complex<double>>;
} // namespace pintr
