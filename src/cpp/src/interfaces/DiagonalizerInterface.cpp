#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <spdlog/spdlog.h>

namespace pairinteraction {
template <typename Scalar>
DiagonalizerInterface<Scalar>::DiagonalizerInterface(FPP fpp) : fpp(fpp) {}

template <typename Scalar>
template <typename ScalarLim>
Eigen::MatrixX<ScalarLim>
DiagonalizerInterface<Scalar>::subtract_mean(const Eigen::MatrixX<Scalar> &matrix, real_t &shift,
                                             int precision) const {
    using real_lim_t = typename traits::NumTraits<ScalarLim>::real_t;
    int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    shift = matrix.diagonal().real().mean();
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        (matrix - shift * Eigen::MatrixX<Scalar>::Identity(dim, dim)).template cast<ScalarLim>();

    // Check if the precision is reachable
    real_lim_t max_entry = shifted_matrix.array().abs().maxCoeff();
    int precision_shift = std::floor(
        -std::log10(shift - std::nextafter(shift, std::numeric_limits<real_t>::lowest())));
    int precision_lim = std::floor(-std::log10(
        max_entry - std::nextafter(max_entry, std::numeric_limits<real_lim_t>::lowest())));
    if (precision > std::min(precision_shift, precision_lim)) {
        SPDLOG_WARN("Because the floating point precision is too low, the energies cannot be "
                    "calculated with a precision of 1e-{} Hartree.",
                    precision);
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
                                    real_t min_eigenvalue, real_t max_eigenvalue,
                                    int precision) const {
    int dim = matrix.rows();

    auto eigensys = eigh(matrix, precision);

    auto *it_begin = std::lower_bound(eigensys.eigenvalues.data(),
                                      eigensys.eigenvalues.data() + dim, min_eigenvalue);
    auto *it_end = std::upper_bound(eigensys.eigenvalues.data(), eigensys.eigenvalues.data() + dim,
                                    max_eigenvalue);
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
                                             int precision) const;
template Eigen::MatrixX<std::complex<float>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, int precision) const;

template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<float> &shifted_eigenvalues,
                                        double shift) const;
template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<float> &shifted_eigenvalues, double shift) const;

template Eigen::MatrixX<double>
DiagonalizerInterface<double>::subtract_mean(const Eigen::MatrixX<double> &matrix, double &shift,
                                             int precision) const;
template Eigen::MatrixX<std::complex<double>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, int precision) const;

template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<double> &shifted_eigenvalues,
                                        double shift) const;
template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<double> &shifted_eigenvalues, double shift) const;
} // namespace pairinteraction
