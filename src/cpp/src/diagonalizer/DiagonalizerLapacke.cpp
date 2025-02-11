#include "pairinteraction/diagonalizer/DiagonalizerLapacke.hpp"

#include "pairinteraction/enums/FPP.hpp"
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

template <typename Scalar, typename ScalarRstr>
EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                   int precision) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using real_rstr_t = typename traits::NumTraits<ScalarRstr>::real_t;
    int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    const real_t shift = matrix.diagonal().real().mean();
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> identity(dim, dim);
    identity.setIdentity();
    Eigen::MatrixX<ScalarRstr> evecs = (matrix - shift * identity).template cast<ScalarRstr>();

    // Check if the precision is reachable
    real_rstr_t max_entry = evecs.array().abs().maxCoeff();
    int precision_shift = std::floor(
        -std::log10(shift - std::nextafter(shift, std::numeric_limits<real_t>::lowest())));
    int precision_rstr = std::floor(-std::log10(
        max_entry - std::nextafter(max_entry, std::numeric_limits<real_rstr_t>::lowest())));
    if (precision > std::min(precision_shift, precision_rstr)) {
        SPDLOG_WARN("Because the floating point precision is too low, the energies cannot be "
                    "calculated with a precision of 1e-{} Hartree.",
                    precision);
    }

    // Diagonalize the shifted matrix
    lapack_int info{}; // will contain return codes
    char jobz = 'V';   // eigenvalues and eigenvectors are computed
    char uplo = 'U';   // full matrix is stored, upper is used

    Eigen::VectorX<real_rstr_t> evals(dim);
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

    // Add the mean of the diagonal elements back to the eigenvalues
    Eigen::VectorX<real_t> eigenvalues = evals.template cast<real_t>();
    eigenvalues.array() += shift;

    return {evecs.sparseView(1, std::pow(10, -precision)).template cast<Scalar>(), eigenvalues};
}

template <typename Scalar>
DiagonalizerLapacke<Scalar>::DiagonalizerLapacke(FPP fpp) : DiagonalizerInterface<Scalar>(fpp) {}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerLapacke<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                  int precision) const {
    switch (this->fpp) {
    case FPP::FLOAT32:
        return dispatch_eigh<Scalar, traits::restricted_t<Scalar, FPP::FLOAT32>>(matrix, precision);
    case FPP::FLOAT64:
        return dispatch_eigh<Scalar, traits::restricted_t<Scalar, FPP::FLOAT64>>(matrix, precision);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

#else

template <typename Scalar>
DiagonalizerLapacke<Scalar>::DiagonalizerLapacke(FPP fpp) : DiagonalizerInterface<Scalar>(fpp) {
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
template class DiagonalizerLapacke<double>;
template class DiagonalizerLapacke<std::complex<double>>;
} // namespace pairinteraction
