// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerLapackeEvr : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    DiagonalizerLapackeEvr(FloatType float_type = FloatType::FLOAT64);
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              double rtol) const override;
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              std::optional<real_t> min_eigenvalue,
                              std::optional<real_t> max_eigenvalue, double rtol) const override;

private:
    template <typename ScalarLim>
    EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                       std::optional<real_t> min_eigenvalue,
                                       std::optional<real_t> max_eigenvalue, double rtol) const;
};

extern template class DiagonalizerLapackeEvr<double>;
extern template class DiagonalizerLapackeEvr<std::complex<double>>;
} // namespace pairinteraction
