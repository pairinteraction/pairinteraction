// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/diagonalize/DiagonalizerLapackeEvd.hpp"

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#ifdef WITH_MKL
#include <mkl.h>
#elif WITH_LAPACKE
#include <lapacke.h>
#endif

namespace pairinteraction {
#if defined(WITH_MKL) || defined(WITH_LAPACKE)
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
template <typename ScalarLim>
EigenSystemH<Scalar> DiagonalizerLapackeEvd<Scalar>::dispatch_eigh(
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix, double rtol) const {
    using real_lim_t = typename traits::NumTraits<ScalarLim>::real_t;
    int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    real_t shift{};
    Eigen::MatrixX<ScalarLim> evecs = this->template subtract_mean<ScalarLim>(matrix, shift, rtol);

    // Diagonalize the shifted matrix
    char jobz = 'V'; // eigenvalues and eigenvectors are computed
    char uplo = 'U'; // full matrix is stored, upper is used

    Eigen::VectorX<real_lim_t> evals(dim);
    lapack_int info = evd(LAPACK_COL_MAJOR, jobz, uplo, dim, evecs.data(), dim, evals.data());

    if (info != 0) {
        if (info < 0) {
            throw std::invalid_argument(
                fmt::format("Diagonalization error: Argument {} to the "
                            "lapacke_evd routine had an illegal value (the counting of the "
                            "arguments starts with one). For a documentation of lapacke_evd, see "
                            "https://www.intel.com/content/www/us/en/docs/onemkl/"
                            "developer-reference-c/2025-1/syevd.html.",
                            -info));
        }
        throw std::runtime_error(fmt::format(
            "Diagonalization error: The lapacke_evd routine failed with error code {}.", info));
    }

    return {evecs.sparseView(1, 0.5 * rtol / std::sqrt(dim)).template cast<Scalar>(),
            this->add_mean(evals, shift)};
}

template <typename Scalar>
DiagonalizerLapackeEvd<Scalar>::DiagonalizerLapackeEvd(FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type) {}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerLapackeEvd<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
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

#else

template <typename Scalar>
DiagonalizerLapackeEvd<Scalar>::DiagonalizerLapackeEvd(FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type) {
    throw std::runtime_error("The lapacke_evd routine is not available in this build. Please use a "
                             "different diagonalizer.");
}

template <typename Scalar>
EigenSystemH<Scalar> DiagonalizerLapackeEvd<Scalar>::eigh(
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/, double /*rtol*/) const {
    std::abort(); // can't happen because the constructor throws
}

#endif // WITH_MKL || WITH_LAPACKE

// Explicit instantiations
template class DiagonalizerLapackeEvd<double>;
template class DiagonalizerLapackeEvd<std::complex<double>>;
} // namespace pairinteraction
