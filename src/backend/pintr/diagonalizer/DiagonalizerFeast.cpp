#include "pintr/diagonalizer/DiagonalizerFeast.hpp"

#include "pintr/utils/eigen_assertion.hpp"
#include "pintr/utils/eigen_compat.hpp"

#include <Eigen/Dense>

#ifdef WITH_MKL
#include <fmt/core.h>
#include <mkl.h>
#endif // WITH_MKL

namespace pintr {
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
DiagonalizerFeast<Scalar>::DiagonalizerFeast(int m0) : m0(m0) {
    if (m0 <= 0) {
        throw std::invalid_argument("The size of the initial subspace m0 (i.e., the number of "
                                    "maximally obtainable eigenvalues) must be positive.");
    }
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                real_t min_eigenvalue, real_t max_eigenvalue, int precision) const {
    int dim = matrix.rows();
    int m0 = std::min(dim, this->m0);

    Eigen::MatrixX<Scalar> hamiltonian = matrix;
    Eigen::VectorX<real_t> evals(dim);
    Eigen::MatrixX<Scalar> evecs(dim, m0); // the first m columns will contain the eigenvectors

    std::vector<MKL_INT> fpm(128);
    feastinit(fpm.data());
    fpm[0] = 0;                  // disable terminal output
    fpm[1] = 5;                  // number of contour points
    fpm[26] = 0;                 // disables matrix checker
    fpm[3] = precision;          // precision
    fpm[5] = 0;                  // user initial subspace
    MKL_INT m{};                 // will contain the number of eigenvalues
    std::vector<real_t> e(m0);   // will contain the first m eigenvalues
    char uplo = 'F';             // full matrix is stored
    MKL_INT info{};              // will contain return codes
    real_t epsout{};             // will contain relative error
    MKL_INT loop{};              // will contain number of used refinement
    std::vector<real_t> res(m0); // will contain the residual errors

    feast(&uplo, &dim, hamiltonian.data(), &dim, fpm.data(), &epsout, &loop, &min_eigenvalue,
          &max_eigenvalue, &m0, evals.data(), evecs.data(), &m, res.data(), &info);

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
                fmt::format("Diagonalization error: The {}-th argument to the FEAST *interface* "
                            "had an illegal value (counting starts at one).",
                            -info - 100));
        }
        if (info >= 100) {
            throw std::invalid_argument(
                fmt::format("Diagonalization error: The {}-th argument to the FEAST "
                            "*initialization* had an illegal value (counting starts at one).",
                            info - 100));
        }
        throw std::runtime_error(
            "Diagonalization error: The FEAST routine failed for an unknown reason.");
    }

    // Restrict to the first m eigenvectors and eigenvalues because the rest is not calculated
    evecs.conservativeResize(dim, m);
    evals.conservativeResize(m);

    return {evecs.sparseView(std::pow(10, -precision), 1), evals};
}

#else

template <typename Scalar>
DiagonalizerFeast<Scalar>::DiagonalizerFeast(int m0) : m0(m0) {
    throw std::runtime_error(
        "The FEAST routine is not available in this build. Please use a different diagonalizer.");
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerFeast<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> & /*matrix*/,
                                real_t /*min_eigenvalue*/, real_t /*max_eigenvalue*/,
                                int /*precision*/) const {
    std::abort(); // can't happen because the constructor throws
}

#endif // WITH_MKL

template class DiagonalizerFeast<float>;
template class DiagonalizerFeast<double>;
template class DiagonalizerFeast<std::complex<float>>;
template class DiagonalizerFeast<std::complex<double>>;
} // namespace pintr
