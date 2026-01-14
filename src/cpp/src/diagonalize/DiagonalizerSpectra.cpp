// SPDX-FileCopyrightText: 2026 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/diagonalize/DiagonalizerSpectra.hpp"

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Spectra/HermEigsSolver.h>
#include <Spectra/MatOp/DenseHermMatProd.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <cassert>
#include <optional>
#include <type_traits>

namespace pairinteraction {

template <typename Scalar>
DiagonalizerSpectra<Scalar>::DiagonalizerSpectra(FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type) {}

template <typename Scalar>
template <typename ScalarLim>
EigenSystemH<Scalar> DiagonalizerSpectra<Scalar>::dispatch_eigh(
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix, std::optional<Eigen::Index> nev,
    std::optional<Eigen::Index> ncv, double rtol) const {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    const int dim = matrix.rows();

    // Subtract the mean of the diagonal elements from the diagonal
    real_t shift{};
    Eigen::MatrixX<ScalarLim> shifted_matrix =
        this->template subtract_mean<ScalarLim>(matrix, shift, rtol);

    // default to half the spectrum
    const int half = std::max(dim / 2, 1);
    const int full = std::min(2 * half, dim);

    // Diagonalize the shifted matrix
    using OpType = std::conditional_t<traits::NumTraits<Scalar>::is_complex_v,
                                      Spectra::DenseHermMatProd<ScalarLim>,
                                      Spectra::DenseSymMatProd<ScalarLim>>;
    using SolType =
        std::conditional_t<traits::NumTraits<Scalar>::is_complex_v, Spectra::HermEigsSolver<OpType>,
                           Spectra::SymEigsSolver<OpType>>;
    OpType op(shifted_matrix);
    SolType eigensolver(op, nev.value_or(half), ncv.value_or(full));
    eigensolver.init();
    eigensolver.compute(Spectra::SortRule::SmallestMagn);
    assert(eigensolver.info() == Spectra::CompInfo::Successful);

    return {eigensolver.eigenvectors()
                .sparseView(1, 0.5 * rtol / std::sqrt(dim))
                .template cast<Scalar>(),
            this->add_mean(eigensolver.eigenvalues(), shift)};
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerSpectra<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                  double rtol) const {
    return this->eigh(matrix, std::nullopt, std::nullopt, rtol);
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerSpectra<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                  std::optional<Eigen::Index> nev, std::optional<Eigen::Index> ncv,
                                  double rtol) const {
    switch (this->float_type) {
    case FloatType::FLOAT32:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT32>>(matrix, nev, ncv,
                                                                               rtol);
    case FloatType::FLOAT64:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT64>>(matrix, nev, ncv,
                                                                               rtol);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

// Explicit instantiations
template class DiagonalizerSpectra<double>;
template class DiagonalizerSpectra<std::complex<double>>;
} // namespace pairinteraction
