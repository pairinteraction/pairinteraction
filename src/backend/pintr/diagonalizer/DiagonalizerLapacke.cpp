#include "pintr/diagonalizer/DiagonalizerLapacke.hpp"

#include "pintr/utils/eigen_assertion.hpp"
#include "pintr/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <fmt/core.h>
#include <mkl_lapacke.h>

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
                                  Range<real_t> allowed_range_of_evals, int precision) const {
    int dim = matrix.rows();

    MKL_INT info;    // will contain return codes
    char jobz = 'V'; // eigenvalues and eigenvectors are computed
    char uplo = 'U'; // full matrix is stored, upper is used

    Eigen::VectorX<real_t> evals(dim);
    Eigen::MatrixX<Scalar> evecs = matrix;

    this->smph_cpu_cores.acquire();
    info = evd(LAPACK_COL_MAJOR, jobz, uplo, dim, evecs.data(), dim, evals.data());
    this->smph_cpu_cores.release();

    if (info != 0) {
        if (info < 0) {
            throw std::runtime_error(fmt::format("Diagonalization error: The {}-th argument to the "
                                                 "LAPACKE routine had an illegal value.",
                                                 -info));
        } else {
            throw std::runtime_error(
                "Diagonalization error: The LAPACK routine failed to compute an eigenvalue.");
        }
    }

    if (allowed_range_of_evals.is_finite()) {
        auto evals_min = allowed_range_of_evals.min;
        auto evals_max = allowed_range_of_evals.max;
        auto it_begin = std::lower_bound(evals.data(), evals.data() + dim, evals_min);
        auto it_end = std::upper_bound(evals.data(), evals.data() + dim, evals_max);
        evecs = evecs.block(0, std::distance(evals.data(), it_begin), dim,
                            std::distance(it_begin, it_end));
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
