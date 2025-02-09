#include "pairinteraction/diagonalizer/DiagonalizerFeast.hpp"

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <cmath>
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

template <typename Scalar, typename ScalarRstr>
EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                   typename traits::NumTraits<Scalar>::real_t min_eigenvalue,
                                   typename traits::NumTraits<Scalar>::real_t max_eigenvalue,
                                   int m0, int precision) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using real_rstr_t = typename traits::NumTraits<ScalarRstr>::real_t;
    int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    const real_t shift = matrix.diagonal().real().mean();
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> identity(dim, dim);
    identity.setIdentity();
    Eigen::MatrixX<ScalarRstr> hamiltonian =
        (matrix - shift * identity).template cast<ScalarRstr>();

    // Check if the precision is reachable
    real_rstr_t entry = hamiltonian.array().abs().maxCoeff();
    int precision_shift = std::floor(
        -std::log10(shift - std::nextafter(shift, std::numeric_limits<real_t>::lowest())));
    int precision_rstr = std::floor(
        -std::log10(entry - std::nextafter(entry, std::numeric_limits<real_rstr_t>::lowest())));
    if (precision > std::min(precision_shift, precision_rstr) +
            1) { // +1 because FEAST needs a higher precision
        SPDLOG_WARN("Because the floating point precision is too low, the energies cannot be "
                    "calculated with a precision of 1e-{} Hartree.",
                    precision);
    }

    // Diagonalize the shifted matrix
    m0 = std::min(dim, m0);

    Eigen::VectorX<real_rstr_t> evals(dim);
    Eigen::MatrixX<ScalarRstr> evecs(dim, m0); // the first m columns will contain the eigenvectors

    std::vector<MKL_INT> fpm(128);
    feastinit(fpm.data());
    fpm[0] = 0;                       // disable terminal output
    fpm[1] = 8;                       // number of contour points
    fpm[4] = 0;                       // do not use initial subspace
    fpm[26] = 0;                      // disables matrix checker
    fpm[2] = precision;               // single precision stopping criteria
    fpm[6] = precision;               // double precision stopping criteria
    MKL_INT m{};                      // will contain the number of eigenvalues
    std::vector<real_rstr_t> e(m0);   // will contain the first m eigenvalues
    char uplo = 'F';                  // full matrix is stored
    MKL_INT info{};                   // will contain return codes
    real_rstr_t epsout{};             // will contain relative error
    MKL_INT loop{};                   // will contain number of used refinement
    std::vector<real_rstr_t> res(m0); // will contain the residual errors
    real_rstr_t min_eigenvalue_rstr = min_eigenvalue - shift;
    real_rstr_t max_eigenvalue_rstr = max_eigenvalue - shift;

    feast(&uplo, &dim, hamiltonian.data(), &dim, fpm.data(), &epsout, &loop, &min_eigenvalue_rstr,
          &max_eigenvalue_rstr, &m0, evals.data(), evecs.data(), &m, res.data(), &info);

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
                            "had an illegal value (counting starts at one).",
                            -info - 100));
        }
        if (info >= 100) {
            throw std::invalid_argument(
                fmt::format("Diagonalization error: Argument {} to the FEAST "
                            "*initialization* had an illegal value (counting starts at one).",
                            info - 100));
        }
        throw std::runtime_error(fmt::format(
            "Diagonalization error: The FEAST routine failed with error code {}.", info));
    }

    // Restrict to the first m eigenvectors and eigenvalues because the rest is not calculated
    evecs.conservativeResize(dim, m);
    evals.conservativeResize(m);

    // Add the mean of the diagonal elements back to the eigenvalues
    Eigen::VectorX<real_t> eigenvalues = evals.template cast<real_t>();
    eigenvalues.array() += shift;

    return {evecs.sparseView(std::pow(10, -precision), 1).template cast<Scalar>(), eigenvalues};
}

template <typename Scalar>
DiagonalizerFeast<Scalar>::DiagonalizerFeast(int m0, FPP fpp)
    : DiagonalizerInterface<Scalar>(fpp), m0(m0) {
    if (m0 <= 0) {
        throw std::invalid_argument("The size of the initial subspace m0 (i.e., the number of "
                                    "maximally obtainable eigenvalues) must be positive.");
    }
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                int /*precision*/) const {
    throw std::invalid_argument("The FEAST routine requires a search interval.");
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                real_t min_eigenvalue, real_t max_eigenvalue, int precision) const {
    switch (this->fpp) {
    case FPP::FLOAT32:
        return dispatch_eigh<Scalar, traits::restricted_t<Scalar, FPP::FLOAT32>>(
            matrix, min_eigenvalue, max_eigenvalue, m0, precision);
    case FPP::FLOAT64:
        return dispatch_eigh<Scalar, traits::restricted_t<Scalar, FPP::FLOAT64>>(
            matrix, min_eigenvalue, max_eigenvalue, m0, precision);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

#else

template <typename Scalar>
DiagonalizerFeast<Scalar>::DiagonalizerFeast(int m0, FPP fpp)
    : DiagonalizerInterface<Scalar>(fpp), m0(m0) {
    throw std::runtime_error(
        "The FEAST routine is not available in this build. Please use a different diagonalizer.");
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                int /*precision*/) const {
    std::abort(); // can't happen because the constructor throws
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                real_t /*min_eigenvalue*/, real_t /*max_eigenvalue*/,
                                int /*precision*/) const {
    std::abort(); // can't happen because the constructor throws
}

#endif // WITH_MKL

// Explicit instantiations
template class DiagonalizerFeast<double>;
template class DiagonalizerFeast<std::complex<double>>;
} // namespace pairinteraction
