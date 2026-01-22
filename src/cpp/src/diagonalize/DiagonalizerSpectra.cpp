// SPDX-FileCopyrightText: 2026 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/diagonalize/DiagonalizerSpectra.hpp"

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Spectra/HermEigsBase.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/Util/SelectionRule.h>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>

namespace pairinteraction {

namespace {
using namespace Spectra;
template <typename OpType>
class HermEigsShiftSolver : public HermEigsBase<OpType, IdentityBOp> {
private:
    using Scalar = typename OpType::Scalar;
    using RealScalar = typename Eigen::NumTraits<Scalar>::Real;
    using Index = Eigen::Index;
    using Array = Eigen::Array<Scalar, Eigen::Dynamic, 1>;

    using Base = HermEigsBase<OpType, IdentityBOp>;
    using Base::m_nev;
    using Base::m_ritz_val;

    const RealScalar m_sigma;

    void sort_ritzpair(SortRule sort_rule) override {
        m_ritz_val.head(m_nev).array() = RealScalar(1) / m_ritz_val.head(m_nev).array() + m_sigma;
        Base::sort_ritzpair(sort_rule);
    }

public:
    HermEigsShiftSolver(OpType &op, Index nev, Index ncv, const RealScalar &sigma)
        : Base(op, IdentityBOp(), nev, ncv), m_sigma(sigma) {
        op.set_shift(m_sigma);
    }
};

} // namespace

template <typename Scalar>
DiagonalizerSpectra<Scalar>::DiagonalizerSpectra(FloatType float_type)
    : DiagonalizerInterface<Scalar>(float_type) {}

template <typename Scalar>
template <typename ScalarLim>
EigenSystemH<Scalar> DiagonalizerSpectra<Scalar>::dispatch_eigh(
    const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix, std::optional<Eigen::Index> nev,
    std::optional<Eigen::Index> ncv, std::optional<real_t> sigma, double rtol) const {

    // Subtract the mean of the diagonal elements from the diagonal
    Eigen::SparseMatrix<ScalarLim> shifted_matrix = matrix.template cast<ScalarLim>();

    // default to half the spectrum
    const Eigen::Index dim = matrix.rows();
    const Eigen::Index half = std::max(dim / 2, Eigen::Index{1});
    const Eigen::Index full = std::min(2 * nev.value_or(half), dim);

    // Diagonalize the shifted matrix
    using OpType = Spectra::SparseSymShiftSolve<ScalarLim>;
    using SolType = HermEigsShiftSolver<OpType>;
    OpType op(shifted_matrix);
    SolType eigensolver(op, nev.value_or(half), ncv.value_or(full), sigma.value_or(0.0));
    eigensolver.init();
    eigensolver.compute(/* selection */ Spectra::SortRule::LargestMagn,
                        /* maxit */ 1000,
                        /* tol */ 1e-10,
                        /* sorting */ Spectra::SortRule::SmallestAlge);
    assert(eigensolver.info() == Spectra::CompInfo::Successful);

    return {eigensolver.eigenvectors()
                .sparseView(1, 0.5 * rtol / std::sqrt(dim))
                .template cast<Scalar>(),
            eigensolver.eigenvalues().template cast<real_t>()};
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerSpectra<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                  double rtol) const {
    return this->eigh(matrix, std::nullopt, std::nullopt, std::nullopt, rtol);
}

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerSpectra<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                  std::optional<real_t> min_eigenvalue,
                                  std::optional<real_t> max_eigenvalue, double rtol) const {
    if (!min_eigenvalue.has_value() || !max_eigenvalue.has_value()) {
        throw std::invalid_argument("The Spectra routine requires a search interval.");
    }
    const real_t sigma = .5 * (min_eigenvalue.value() + max_eigenvalue.value());
    const auto [global_lower, global_upper, lower_count, upper_count] =
        this->gershgorin_bounds(matrix, min_eigenvalue.value(), max_eigenvalue.value());

    auto eigensys = this->eigh(matrix, upper_count, std::nullopt, sigma, rtol);

    const int dim = matrix.rows();
    const int nev = eigensys.eigenvalues.rows();

    auto *it_begin = std::lower_bound(eigensys.eigenvalues.data(),
                                      eigensys.eigenvalues.data() + nev, min_eigenvalue.value());
    auto *it_end = std::upper_bound(eigensys.eigenvalues.data(), eigensys.eigenvalues.data() + nev,
                                    max_eigenvalue.value());
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

template <typename Scalar>
EigenSystemH<Scalar>
DiagonalizerSpectra<Scalar>::eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                  std::optional<Eigen::Index> nev, std::optional<Eigen::Index> ncv,
                                  std::optional<real_t> sigma, double rtol) const {
    switch (this->float_type) {
    case FloatType::FLOAT32:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT32>>(matrix, nev, ncv,
                                                                               sigma, rtol);
    case FloatType::FLOAT64:
        return dispatch_eigh<traits::restricted_t<Scalar, FloatType::FLOAT64>>(matrix, nev, ncv,
                                                                               sigma, rtol);
    default:
        throw std::invalid_argument("Unsupported floating point precision.");
    }
}

// Explicit instantiations
template class DiagonalizerSpectra<double>;
template class DiagonalizerSpectra<std::complex<double>>;
} // namespace pairinteraction
