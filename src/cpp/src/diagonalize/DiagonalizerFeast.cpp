// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/diagonalize/DiagonalizerFeast.hpp"

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <optional>
#include <spdlog/spdlog.h>

#ifdef WITH_MKL
#include <fmt/core.h>
#include <mkl.h>
#endif // WITH_MKL

namespace pairinteraction {
#ifdef WITH_MKL
void feast(const char *uplo, const MKL_INT *n, const float *a, const MKL_INT *lda, MKL_INT *fpm,
           float *epsout, MKL_INT *loop, const float *emin, const float *emax, MKL_INT *m0,
           float *e, float *x, MKL_INT *m, float *res, MKL_INT *info) {
    sfeast_syev(uplo, n, a, lda, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info);
}

void feast(const char *uplo, const MKL_INT *n, const double *a, const MKL_INT *lda, MKL_INT *fpm,
           double *epsout, MKL_INT *loop, const double *emin, const double *emax, MKL_INT *m0,
           double *e, double *x, MKL_INT *m, double *res, MKL_INT *info) {
    dfeast_syev(uplo, n, a, lda, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info);
}

void feast(const char *uplo, const MKL_INT *n, const MKL_Complex8 *a, const MKL_INT *lda,
           MKL_INT *fpm, float *epsout, MKL_INT *loop, const float *emin, const float *emax,
           MKL_INT *m0, float *e, MKL_Complex8 *x, MKL_INT *m, float *res, MKL_INT *info) {
    cfeast_heev(uplo, n, a, lda, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info);
}

void feast(const char *uplo, const MKL_INT *n, const MKL_Complex16 *a, const MKL_INT *lda,
           MKL_INT *fpm, double *epsout, MKL_INT *loop, const double *emin, const double *emax,
           MKL_INT *m0, double *e, MKL_Complex16 *x, MKL_INT *m, double *res, MKL_INT *info) {
    zfeast_heev(uplo, n, a, lda, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info);
}

template <typename Scalar>
template <typename ScalarLim>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                         real_t min_eigenvalue, real_t max_eigenvalue,
                                         double rtol) const {
    using real_lim_t = typename traits::NumTraits<ScalarLim>::real_t;
    int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    real_t shift{};
    Eigen::MatrixX<ScalarLim> hamiltonian =
        this->template subtract_mean<ScalarLim>(matrix, shift, rtol);

    // Diagonalize the shifted matrix
    int m0 = std::min(dim, this->m0);

    Eigen::VectorX<real_lim_t> evals(dim);
    Eigen::MatrixX<ScalarLim> evecs(dim, m0); // the first m columns will contain the eigenvectors

    double targeted_trace_relative_error = rtol;
    int precision_feast = static_cast<int>(std::ceil(-std::log10(targeted_trace_relative_error)));

    std::vector<MKL_INT> fpm(128);
    feastinit(fpm.data());
    fpm[0] = 0;                      // disable terminal output
    fpm[1] = 8;                      // number of contour points
    fpm[4] = 0;                      // do not use initial subspace
    fpm[26] = 0;                     // disables matrix checker
    fpm[2] = precision_feast;        // single precision stopping criteria
    fpm[6] = precision_feast;        // double precision stopping criteria
    MKL_INT m{};                     // will contain the number of eigenvalues
    std::vector<real_lim_t> e(m0);   // will contain the first m eigenvalues
    char uplo = 'F';                 // full matrix is stored
    MKL_INT info{};                  // will contain return codes
    real_lim_t epsout{};             // will contain relative error
    MKL_INT loop{};                  // will contain number of used refinement
    std::vector<real_lim_t> res(m0); // will contain the residual errors
    real_lim_t min_eigenvalue_lim = min_eigenvalue - shift;
    real_lim_t max_eigenvalue_lim = max_eigenvalue - shift;

    feast(&uplo, &dim, hamiltonian.data(), &dim, fpm.data(), &epsout, &loop, &min_eigenvalue_lim,
          &max_eigenvalue_lim, &m0, evals.data(), evecs.data(), &m, res.data(), &info);

    // https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/extended-eigensolver-output-details.html
    if (info != 0) {
        if (info == 202) {
            throw std::invalid_argument(
                "Diagonalization error: Problem with size of the system n (n≤0).");
        }
        if (info == 201) {
            throw std::invalid_argument(
                "Diagonalization error: Problem with size of initial subspace m0 (m0≤0 or m0>n).");
        }
        if (info == 200) {
            throw std::invalid_argument(
                "Diagonalization error: Problem with emin,emax (emin≥emax).");
        }
        if (info == 3) {
            throw std::invalid_argument(
                "Diagonalization error: Size of the subspace m0 is too small (m0<m).");
        }
        if (info == 2) {
            throw std::runtime_error(
                "Diagonalization error: No convergence (number of iteration loops >fpm[3]).");
        }
        if (info == 1) {
            throw std::runtime_error(
                "Diagonalization error: No eigenvalue found in the search interval.");
        }
        if (info == -1) {
            throw std::runtime_error(
                "Diagonalization error: Internal error for allocation memory.");
        }
        if (info == -2) {
            throw std::runtime_error(
                "Diagonalization error: Internal error of the inner system solver. Possible "
                "reasons: not enough memory for inner linear system solver or inconsistent input.");
        }
        if (info == -3) {
            throw std::runtime_error(
                "Diagonalization error: Internal error of the reduced eigenvalue solver. Possible "
                "cause: matrix may not be positive definite.");
        }
        if (info == -4) {
            throw std::invalid_argument("Diagonalization error: Matrix is not positive definite.");
        }
        if (info <= 100) {
            throw std::invalid_argument(
                fmt::format("Diagonalization error: Argument {} to the FEAST *interface* "
                            "had an illegal value (the counting of the arguments starts with one). "
                            "For a documentation of the feast interface, see "
                            "https://www.intel.com/content/www/us/en/docs/onemkl/"
                            "developer-reference-c/2025-1/feast-syev-feast-heev.html.",
                            -info - 100));
        }
        if (info >= 100) {
            throw std::invalid_argument(fmt::format(
                "Diagonalization error: Argument {} to the FEAST "
                "*initialization* had an illegal value (the counting of the arguments starts with "
                "one). For a documentation of the feast initialization, see "
                "https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2025-1/"
                "extended-eigensolver-input-parameters.html.",
                info - 100));
        }
        throw std::runtime_error(fmt::format(
            "Diagonalization error: The FEAST routine failed with error code {}.", info));
    }

    // Restrict to the first m eigenvectors and eigenvalues because the rest is not calculated
    evecs.conservativeResize(dim, m);
    evals.conservativeResize(m);

    return {evecs.sparseView(1, 0.5 * rtol / std::sqrt(dim)).template cast<Scalar>(),
            this->add_mean(evals, shift)};
}

template <typename Scalar>
DiagonalizerFeast<Scalar>::DiagonalizerFeast(int m0, FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type), m0(m0) {
    if (m0 <= 0) {
        throw std::invalid_argument("The size of the initial subspace m0 (i.e., the number of "
                                    "maximally obtainable eigenvalues) must be positive.");
    }
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                double /*rtol*/) const {
    throw std::invalid_argument("The FEAST routine requires a search interval.");
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                std::optional<real_t> min_eigenvalue,
                                std::optional<real_t> max_eigenvalue, double rtol) const {
    if (!min_eigenvalue.has_value() || !max_eigenvalue.has_value()) {
        throw std::invalid_argument("The FEAST routine requires a search interval.");
    }
    switch (this->float_type) {
    case FloatType::FLOAT32:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT32>>(
            matrix, min_eigenvalue.value(), max_eigenvalue.value(), rtol);
    case FloatType::FLOAT64:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT64>>(
            matrix, min_eigenvalue.value(), max_eigenvalue.value(), rtol);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

#else

template <typename Scalar>
DiagonalizerFeast<Scalar>::DiagonalizerFeast(int m0, FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type), m0(m0) {
    throw std::runtime_error(
        "The FEAST routine is not available in this build. Please use a different diagonalizer.");
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                double /*rtol*/) const {
    std::abort(); // can't happen because the constructor throws
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                std::optional<real_t> /*min_eigenvalue*/,
                                std::optional<real_t> /*max_eigenvalue*/, double /*rtol*/) const {
    std::abort(); // can't happen because the constructor throws
}

#endif // WITH_MKL

// Explicit instantiations
template class DiagonalizerFeast<double>;
template class DiagonalizerFeast<std::complex<double>>;
} // namespace pairinteraction
