#include "pintr/diagonalizer/DiagonalizerFeast.hpp"

#include "pintr/utils/eigen_assertion.hpp"
#include "pintr/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <limits>

namespace pintr {
template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                real_t min_eigenvalue, real_t max_eigenvalue, int precision) const {
    int dim = matrix.rows();

    Eigen::VectorX<real_t> evals(dim);
    Eigen::MatrixX<Scalar> evecs = matrix;

    // TODO
    (void)min_eigenvalue;
    (void)max_eigenvalue;

    return {evecs.sparseView(std::pow(10, -precision), 1), evals};
}

template class DiagonalizerFeast<float>;
template class DiagonalizerFeast<double>;
template class DiagonalizerFeast<std::complex<float>>;
template class DiagonalizerFeast<std::complex<double>>;
} // namespace pintr
