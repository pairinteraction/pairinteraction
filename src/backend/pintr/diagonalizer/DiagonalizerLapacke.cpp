#include "pintr/diagonalizer/DiagonalizerLapacke.hpp"

#include "pintr/utils/eigen_assertion.hpp"
#include "pintr/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <fmt/core.h>

#ifdef WITH_MKL
#include <mkl.h>
#else
#include <lapacke.h>
#endif

namespace pintr {
lapack_int evd(int matrix_layout, char jobz, char uplo, lapack_int n, float *a, lapack_int lda,
               float *w) {
    return LAPACKE_ssyevd(matrix_layout, jobz, uplo, n, a, lda, w);
};

lapack_int evd(int matrix_layout, char jobz, char uplo, lapack_int n, double *a, lapack_int lda,
               double *w) {
    return LAPACKE_dsyevd(matrix_layout, jobz, uplo, n, a, lda, w);
};

lapack_int evd(int matrix_layout, char jobz, char uplo, lapack_int n, lapack_complex_float *a,
               lapack_int lda, float *w) {
    return LAPACKE_cheevd(matrix_layout, jobz, uplo, n, a, lda, w);
};

lapack_int evd(int matrix_layout, char jobz, char uplo, lapack_int n, lapack_complex_double *a,
               lapack_int lda, double *w) {
    return LAPACKE_zheevd(matrix_layout, jobz, uplo, n, a, lda, w);
};

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerLapacke<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                  real_t min_eigenvalue, real_t max_eigenvalue,
                                  int precision) const {
    int dim = matrix.rows();

    lapack_int info{}; // will contain return codes
    char jobz = 'V';   // eigenvalues and eigenvectors are computed
    char uplo = 'U';   // full matrix is stored, upper is used

    Eigen::VectorX<real_t> evals(dim);
    Eigen::MatrixX<Scalar> evecs = matrix;

    info = evd(LAPACK_COL_MAJOR, jobz, uplo, dim, evecs.data(), dim, evals.data());

    if (info != 0) {
        if (info < 0) {
            throw std::invalid_argument(
                fmt::format("Diagonalization error: The {}-th argument to the "
                            "LAPACKE routine had an illegal value.",
                            -info));
        }
        throw std::runtime_error(
            "Diagonalization error: The LAPACK routine failed to compute an eigenvalue.");
    }

    if (min_eigenvalue > std::numeric_limits<real_t>::lowest() / 2 ||
        max_eigenvalue < std::numeric_limits<real_t>::max() / 2) {
        auto *it_begin = std::lower_bound(evals.data(), evals.data() + dim, min_eigenvalue);
        auto *it_end = std::upper_bound(evals.data(), evals.data() + dim, max_eigenvalue);
        evecs = evecs
                    .block(0, std::distance(evals.data(), it_begin), dim,
                           std::distance(it_begin, it_end))
                    .eval();
        evals =
            evals.segment(std::distance(evals.data(), it_begin), std::distance(it_begin, it_end));
    }

    return {evecs.sparseView(std::pow(10, -precision), 1), evals};
}

template class DiagonalizerLapacke<float>;
template class DiagonalizerLapacke<double>;
template class DiagonalizerLapacke<std::complex<float>>;
template class DiagonalizerLapacke<std::complex<double>>;
} // namespace pintr
