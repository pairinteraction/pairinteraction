#include "pairinteraction/diagonalizer/DiagonalizerLapacke.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <fmt/core.h>

#ifdef WITH_MKL
#include <mkl.h>
#elif WITH_LAPACKE
#include <lapacke.h>
#endif

namespace pairinteraction {
#if defined WITH_MKL || defined WITH_LAPACKE
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
DiagonalizerLapacke<Scalar>::DiagonalizerLapacke() = default;

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerLapacke<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
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
            throw std::invalid_argument(fmt::format("Diagonalization error: Argument {} to the "
                                                    "LAPACKE routine had an illegal value.",
                                                    -info));
        }
        throw std::runtime_error(fmt::format(
            "Diagonalization error: The LAPACK routine failed with error code {}.", info));
    }

    return {evecs.sparseView(std::pow(10, -precision), 1), evals};
}

#else

template <typename Scalar>
DiagonalizerLapacke<Scalar>::DiagonalizerLapacke() {
    throw std::runtime_error(
        "The LAPACKE routine is not available in this build. Please use a different diagonalizer.");
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerLapacke<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                  int /*precision*/) const {
    std::abort(); // can't happen because the constructor throws
}

#endif // WITH_MKL || WITH_LAPACKE

// Explicit instantiations
template class DiagonalizerLapacke<float>;
template class DiagonalizerLapacke<double>;
template class DiagonalizerLapacke<std::complex<float>>;
template class DiagonalizerLapacke<std::complex<double>>;
} // namespace pairinteraction
