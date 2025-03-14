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
                                             double atol) const {
    using real_lim_t = typename traits::NumTraits<ScalarLim>::real_t;
    int dim = matrix.rows();

    assert(atol > 0);
    assert(atol < 1);

    // Subtract the mean of the diagonal elements from the diagonal
    shift = matrix.diagonal().real().mean();
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        (matrix - shift * Eigen::MatrixX<Scalar>::Identity(dim, dim)).template cast<ScalarLim>();

    ///////////////////////////////////////
    // Ensure accuracy of the eigenvectors
    ///////////////////////////////////////

    // We are interested in the maximum error of an entry in an eigenvector that can result from the
    // limited floating point precision. Thus, we calculate the error of the largest possible entry,
    // which is 1 (the eigenvectors are normalized):
    // |1 - nextafter(1)| = eps/2

    // By comparing atol with eps/2, we can estimate whether the floating point error or the error
    // from the specified tolerance dominates.

    double error_from_float_type =
        1e1 * std::numeric_limits<real_lim_t>::epsilon() / 2; // 1e1 is a safety factor

    if (error_from_float_type > atol) {
        SPDLOG_WARN(
            "Because the floating point precision is too low, the "
            "eigenvectors cannot be calculated accurately. The estimated floating point error ({}) "
            "is larger than the specified tolerance ({}).",
            error_from_float_type, atol);
    }

    ///////////////////////////////////////
    // Ensure accuracy of the eigenvalues
    ///////////////////////////////////////

    // Estimate the error of the eigenenergies that would result from calculating the eigenenergies
    // by transforming the Hamiltonian with the eigenvector matrix, assuming that all entries of the
    // eigenvector matrix that are smaller than the tolerance are set to zero:
    // |Tr(D_exact) - Tr(D)|/dim(D) <= 2*atol*||shifted_matrix||

    // On the other hand, we can estimate the error of the eigenenergies that results from the
    // limited floating point precision. As the diagonalization algorithms are typically "backward
    // stable", the computed eigenvalues are typically the eigenvalues of a perturbed matrix
    // "shifted_matrix + E", with ||E|| <= c*eps*||shifted_matrix||, where c is a small constant and
    // eps is the machine epsilon. The error of the eigenvalues is then bounded by:
    // |Tr(D_exact) - Tr(D)|/dim(D) <= c*eps*||shifted_matrix||.

    // Thus, by comparing 2*atol with c*eps, we can estimate whether the floating point error or the
    // error from the eigenvectors dominates. Notably, for c=10, this comparison is equivalent to
    // the check that we performed above! Thus, we can skip this comparison.

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
                                    std::optional<real_t> max_eigenvalue, double atol) const {
    int dim = matrix.rows();

    auto eigensys = eigh(matrix, atol);

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
                                             double atol) const;
template Eigen::MatrixX<std::complex<float>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, double atol) const;

template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<float> &shifted_eigenvalues,
                                        double shift) const;
template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<float> &shifted_eigenvalues, double shift) const;

template Eigen::MatrixX<double>
DiagonalizerInterface<double>::subtract_mean(const Eigen::MatrixX<double> &matrix, double &shift,
                                             double atol) const;
template Eigen::MatrixX<std::complex<double>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, double atol) const;

template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<double> &shifted_eigenvalues,
                                        double shift) const;
template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<double> &shifted_eigenvalues, double shift) const;
} // namespace pairinteraction
