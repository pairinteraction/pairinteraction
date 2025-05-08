// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/diagonalize/DiagonalizerLapackeEvr.hpp"

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <fmt/core.h>
#include <limits>
#include <spdlog/spdlog.h>

#ifdef WITH_MKL
#include <mkl.h>
#elif WITH_LAPACKE
#include <lapacke.h>
#endif

namespace pairinteraction {
#if defined(WITH_MKL) || defined(WITH_LAPACKE)
lapack_int evr(int matrix_layout, char jobz, char range, char uplo, lapack_int n, float *a,
               lapack_int lda, float vl, float vu, lapack_int il, lapack_int iu, float abstol,
               lapack_int *m, float *w, float *z, lapack_int ldz, lapack_int *isuppz) {
    return LAPACKE_ssyevr(matrix_layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,
                          z, ldz, isuppz);
};

lapack_int evr(int matrix_layout, char jobz, char range, char uplo, lapack_int n, double *a,
               lapack_int lda, double vl, double vu, lapack_int il, lapack_int iu, double abstol,
               lapack_int *m, double *w, double *z, lapack_int ldz, lapack_int *isuppz) {
    return LAPACKE_dsyevr(matrix_layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,
                          z, ldz, isuppz);
};

lapack_int evr(int matrix_layout, char jobz, char range, char uplo, lapack_int n,
               lapack_complex_float *a, lapack_int lda, float vl, float vu, lapack_int il,
               lapack_int iu, float abstol, lapack_int *m, float *w, lapack_complex_float *z,
               lapack_int ldz, lapack_int *isuppz) {
    return LAPACKE_cheevr(matrix_layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,
                          z, ldz, isuppz);
};

lapack_int evr(int matrix_layout, char jobz, char range, char uplo, lapack_int n,
               lapack_complex_double *a, lapack_int lda, double vl, double vu, lapack_int il,
               lapack_int iu, double abstol, lapack_int *m, double *w, lapack_complex_double *z,
               lapack_int ldz, lapack_int *isuppz) {
    return LAPACKE_zheevr(matrix_layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w,
                          z, ldz, isuppz);
};

template <typename Scalar>
template <typename ScalarLim>
EigenSystemH<Scalar> DiagonalizerLapackeEvr<Scalar>::dispatch_eigh(
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
    std::optional<real_t> min_eigenvalue, std::optional<real_t> max_eigenvalue, double rtol) const {
    using real_lim_t = typename traits::NumTraits<ScalarLim>::real_t;
    int dim = matrix.rows();
    bool has_range = min_eigenvalue.has_value() || max_eigenvalue.has_value();

    // Subtract the mean of the diagonal elements from the diagonal
    real_t shift{};
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        this->template subtract_mean<ScalarLim>(matrix, shift, rtol);
    real_lim_t scaling = shifted_matrix.norm();

    if (scaling == 0) {
        scaling = 1; // Avoid division by zero if the matrix is zero
    }

    shifted_matrix /= scaling; // This seems to increase the numerical stability of lapacke_evr

    // Diagonalize the shifted matrix
    lapack_int m = 0;                        // Number of eigenvalues found
    char jobz = 'V';                         // Compute eigenvectors
    char range_char = has_range ? 'V' : 'A'; // Compute all eigenvalues if no range is specified
    char uplo = 'U';                         // Matrix is stored in upper-triangular part
    real_lim_t abstol = 0.1 * rtol; // 0.1 is a safety factor, abstol ~ rtol because ||H||=1
    lapack_int il = 0;              // Lower index bound if 'I'
    lapack_int iu = 0;              // Upper index bound if 'I'
    real_lim_t vl = (min_eigenvalue.value_or(-1) - shift) / scaling; // Lower eval bounds if 'V'
    real_lim_t vu = (max_eigenvalue.value_or(1) - shift) / scaling;  // Upper eval bounds if 'V'

    Eigen::VectorX<real_lim_t> evals(dim);                        // Eigenvalues
    Eigen::MatrixX<ScalarLim> evecs(dim, dim);                    // Eigenvectors
    std::vector<lapack_int> isuppz(static_cast<size_t>(2 * dim)); // Workspace
    lapack_int info =
        evr(LAPACK_COL_MAJOR, jobz, range_char, uplo, dim, shifted_matrix.data(), dim, vl, vu, il,
            iu, abstol, &m, evals.data(), evecs.data(), dim, isuppz.data());

    if (info != 0) {
        if (info < 0) {
            throw std::invalid_argument(
                fmt::format("Diagonalization error: Argument {} to the "
                            "lapacke_evr routine had an illegal value (the counting of the "
                            "arguments starts with one). For a documentation of lapacke_evr, see "
                            "https://www.intel.com/content/www/us/en/docs/onemkl/"
                            "developer-reference-c/2025-1/syevr.html.",
                            -info));
        }
        throw std::runtime_error(
            fmt::format("Diagonalization error: The lapacke_evr routine failed with error code {}. "
                        "Try to set a lower 'rtol'.",
                        info));
    }

    // Restrict to the first m eigenvectors and eigenvalues because the rest is not calculated
    evals.conservativeResize(m);
    evals *= scaling;

    return {evecs.leftCols(m).sparseView(1, 0.5 * rtol / std::sqrt(dim)).template cast<Scalar>(),
            this->add_mean(evals, shift)};
}

template <typename Scalar>
DiagonalizerLapackeEvr<Scalar>::DiagonalizerLapackeEvr(FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type) {}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerLapackeEvr<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                     double rtol) const {
    return this->eigh(matrix, std::nullopt, std::nullopt, rtol);
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerLapackeEvr<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                     std::optional<real_t> min_eigenvalue,
                                     std::optional<real_t> max_eigenvalue, double rtol) const {
    switch (this->float_type) {
    case FloatType::FLOAT32:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT32>>(
            matrix, min_eigenvalue, max_eigenvalue, rtol);
    case FloatType::FLOAT64:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT64>>(
            matrix, min_eigenvalue, max_eigenvalue, rtol);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

#else

template <typename Scalar>
DiagonalizerLapackeEvr<Scalar>::DiagonalizerLapackeEvr(FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type) {
    throw std::runtime_error("The lapacke_evr routine is not available in this build. Please use a "
                             "different diagonalizer.");
}

template <typename Scalar>
EigenSystemH<Scalar> DiagonalizerLapackeEvr<Scalar>::eigh(
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/, double /*rtol*/) const {
    std::abort(); // can't happen because the constructor throws
}

template <typename Scalar>
EigenSystemH<Scalar> DiagonalizerLapackeEvr<Scalar>::eigh(
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
    std::optional<real_t> /*min_eigenvalue*/, std::optional<real_t> /*max_eigenvalue*/,
    double /*rtol*/) const {
    std::abort(); // can't happen because the constructor throws
}

#endif // WITH_MKL || WITH_LAPACKE

// Explicit instantiations
template class DiagonalizerLapackeEvr<double>;
template class DiagonalizerLapackeEvr<std::complex<double>>;
} // namespace pairinteraction
