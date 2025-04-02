// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <optional>
#include <spdlog/spdlog.h>

namespace pairinteraction {
template <typename Scalar>
DiagonalizerInterface<Scalar>::DiagonalizerInterface(FloatType float_type)
    : float_type(float_type) {}

template <typename Scalar>
template <typename ScalarLim>
Eigen::MatrixX<ScalarLim>
DiagonalizerInterface<Scalar>::subtract_mean(const Eigen::MatrixX<Scalar> &matrix, real_t &shift,
                                             double rtol) const {
    using real_lim_t = typename traits::NumTraits<ScalarLim>::real_t;
    int dim = matrix.rows();

    assert(rtol > 0);
    assert(rtol < 1);

    // Subtract the mean of the diagonal elements from the diagonal
    shift = matrix.diagonal().real().mean();
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        (matrix - shift * Eigen::MatrixX<Scalar>::Identity(dim, dim)).template cast<ScalarLim>();

    // Ensure the accuracy the results
    double floating_point_error = 5 * std::numeric_limits<real_lim_t>::epsilon();

    if (floating_point_error > rtol) {
        SPDLOG_WARN(
            "Because the floating point precision is too low, the eigenvalues cannot be calculated "
            "accurately. The estimated floating point error ({} * ||H||) is larger than the "
            "specified tolerance ({} * ||H||). Try to use a 'float_type' with higher precision or "
            "a larger 'rtol'.",
            floating_point_error, rtol);
    }

    // Note that this check is not as handwavy as it seems at first glance:
    //
    // Let us first analyze the error of the eigenvalues due to the limited floating point
    // precision. Under the assumption that the diagonalization routine is backward stable and the
    // we know that the calculated eigenvalues and eigenvectors correspond to the diagonalization of
    // a matrix H+E, where H is the original matrix and E a small perturbation with
    // ||E|| <= c*eps*||H||. Here, for an hermitian matrix, c is a modest constant and eps is the
    // machine epsilon. Taking this result, we can use the Bauer-Fike theorem, simplified for an
    // hermitian matrix H, to bound the error of the eigenvalues:
    //
    //     |lambda - lambda_exact| <= c*eps*||H||,
    //
    // where lambda is the calculated eigenvalue and lambda_exact is the closest exact eigenvalue of
    // H. Thus, to ensure that the relative tolerance rtol can be met, we need c*eps*||H|| <=
    // rtol*||H||, which we check in the code above.
    //
    // If we know that for a user, it is permissible that |lambda - lambda_exact| ~ rtol*||H||, we
    // assume that the eigenvectors are allowed to have an error as well. To obtain a reasonable
    // bound for the strength of this error, we demand that the error of an eigenvalue, which would
    // result from calculating the eigenvalue by transforming the Hamiltonian with the
    // eigenvector v, should be smaller than rtol * ||H||. For v = v_exact + w, we obtain:
    //
    //     |v^dag * H * v - lambda_exact| <= 2*||w||*||H||
    //
    // Thus, to keep this error below rtol * ||H||, we need ||w|| <= rtol/2.
    //
    // In the diagonalizer classes derived from the DiagonalizerInterface, we will make use of the
    // fact that a small error of the eigenvectors is tolerable to increase the sparsity of the
    // eigenvector matrix by setting every entry smaller than rtol/(2*sqrt(dim(H))) to zero.
    // The factor 1/sqrt(dim(H)) ensures that the bound is met also if multiple entries are
    // set to zero.

    return shifted_matrix;
}

template <typename Scalar>
template <typename RealLim>
Eigen::VectorX<typename DiagonalizerInterface<Scalar>::real_t>
DiagonalizerInterface<Scalar>::add_mean(const Eigen::VectorX<RealLim> &shifted_eigenvalues,
                                        real_t shift) const {
    // Add the mean of the diagonal elements back to the eigenvalues
    Eigen::VectorX<real_t> eigenvalues = shifted_eigenvalues.template cast<real_t>();
    eigenvalues.array() += shift;
    return eigenvalues;
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerInterface<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                    std::optional<real_t> min_eigenvalue,
                                    std::optional<real_t> max_eigenvalue, double rtol) const {
    int dim = matrix.rows();

    auto eigensys = eigh(matrix, rtol);

    auto *it_begin =
        std::lower_bound(eigensys.eigenvalues.data(), eigensys.eigenvalues.data() + dim,
                         min_eigenvalue.value_or(std::numeric_limits<real_t>::lowest() / 2));
    auto *it_end =
        std::upper_bound(eigensys.eigenvalues.data(), eigensys.eigenvalues.data() + dim,
                         max_eigenvalue.value_or(std::numeric_limits<real_t>::max() / 2));
    eigensys.eigenvectors = eigensys.eigenvectors
                                .block(0, std::distance(eigensys.eigenvalues.data(), it_begin), dim,
                                       std::distance(it_begin, it_end))
                                .eval();
    eigensys.eigenvalues = eigensys.eigenvalues
                               .segment(std::distance(eigensys.eigenvalues.data(), it_begin),
                                        std::distance(it_begin, it_end))
                               .eval();

    return eigensys;
}

// Explicit instantiations
template class DiagonalizerInterface<double>;
template class DiagonalizerInterface<std::complex<double>>;

template Eigen::MatrixX<float>
DiagonalizerInterface<double>::subtract_mean(const Eigen::MatrixX<double> &matrix, double &shift,
                                             double rtol) const;
template Eigen::MatrixX<std::complex<float>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, double rtol) const;

template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<float> &shifted_eigenvalues,
                                        double shift) const;
template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<float> &shifted_eigenvalues, double shift) const;

template Eigen::MatrixX<double>
DiagonalizerInterface<double>::subtract_mean(const Eigen::MatrixX<double> &matrix, double &shift,
                                             double rtol) const;
template Eigen::MatrixX<std::complex<double>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, double rtol) const;

template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<double> &shifted_eigenvalues,
                                        double shift) const;
template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<double> &shifted_eigenvalues, double shift) const;
} // namespace pairinteraction
