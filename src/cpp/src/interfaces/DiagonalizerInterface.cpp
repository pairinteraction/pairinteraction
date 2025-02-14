#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <optional>
#include <spdlog/spdlog.h>

namespace pairinteraction {
template <typename Scalar>
DiagonalizerInterface<Scalar>::DiagonalizerInterface(FPP fpp) : fpp(fpp) {}

template <typename Scalar>
template <typename ScalarLim>
Eigen::MatrixX<ScalarLim>
DiagonalizerInterface<Scalar>::subtract_mean(const Eigen::MatrixX<Scalar> &matrix, real_t &shift,
                                             double atol) const {
    using real_lim_t = typename traits::NumTraits<ScalarLim>::real_t;
    int dim = matrix.rows();

    assert(atol >= 0);
    assert(atol < 1);

    // Subtract the mean of the diagonal elements from the diagonal
    shift = matrix.diagonal().real().mean();
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        (matrix - shift * Eigen::MatrixX<Scalar>::Identity(dim, dim)).template cast<ScalarLim>();

    ///////////////////////////////////////
    // Ensure accuracy of the eigenvectors
    ///////////////////////////////////////

    double floating_point_precision =
        1 - std::nextafter(static_cast<real_lim_t>(1), std::numeric_limits<real_lim_t>::lowest());

    if (floating_point_precision * 1e1 > atol) { // 1e1 is a safety factor
        SPDLOG_WARN(
            "Because the floating point precision is lacking, the "
            "eigenvectors cannot be calculated accurately. The floating point precision ({}) "
            "is similar or worse than the specified tolerance ({}).",
            floating_point_precision, atol);
    }

    ///////////////////////////////////////
    // Accuracy of the eigenvalues
    ///////////////////////////////////////

    // Estimate the error of the eigenenergies due to limited floating point precision
    real_lim_t max_entry = shifted_matrix.array().abs().maxCoeff();
    floating_point_precision =
        max_entry - std::nextafter(max_entry, std::numeric_limits<real_lim_t>::lowest());

    // Estimate the error of the eigenenergies that would result from calculating the eigenenergies
    // by transforming the Hamiltonian with the eigenvector matrix, assuming that all entries of the
    // eigenvector matrix that are smaller than the tolerance are set to zero.
    double precision_from_tolerance =
        2 * atol * shifted_matrix.norm(); // upper bound of |Tr(D_exact) - Tr(D)|/dim(D)

    if (floating_point_precision * 1e1 > precision_from_tolerance) { // 1e1 is a safety factor
        SPDLOG_WARN("Because the floating point precision is lacking, the "
                    "eigenvalues cannot be calculated accurately. The floating point precision ({} "
                    "Hartree) is similar or worse than the precision estimated from the specified "
                    "tolerance ({} Hartree).",
                    floating_point_precision, precision_from_tolerance);
    }

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
