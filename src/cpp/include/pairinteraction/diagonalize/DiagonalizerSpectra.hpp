// SPDX-FileCopyrightText: 2026 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <optional>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerSpectra : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    DiagonalizerSpectra(std::optional<Eigen::Index> ncv = std::nullopt,
                        FloatType float_type = FloatType::FLOAT64);
    EigenSystemH<Scalar> eigh_full(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                   double rtol) const override;
    EigenSystemH<Scalar> eigh_range(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                    std::optional<real_t> min_eigenvalue,
                                    std::optional<real_t> max_eigenvalue,
                                    double rtol) const override;
    EigenSystemH<Scalar>
    eigh_shift_invert(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                      std::optional<Eigen::Index> num_eigenvalues, std::optional<real_t> sigma,
                      double rtol) const override;

private:
    std::optional<Eigen::Index> ncv;
    template <typename ScalarLim>
    EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                       Eigen::Index num_eigenvalues, std::optional<real_t> sigma,
                                       double rtol) const;
};

extern template class DiagonalizerSpectra<double>;
extern template class DiagonalizerSpectra<std::complex<double>>;
} // namespace pairinteraction
