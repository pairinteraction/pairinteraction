#pragma once

#include "pairinteraction/enums/FloatType.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerLapackeEvd : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    DiagonalizerLapackeEvd(FloatType float_type = FloatType::FLOAT64);
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              double atol) const override;

private:
    template <typename ScalarLim>
    EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                       double atol) const;
};

extern template class DiagonalizerLapackeEvd<double>;
extern template class DiagonalizerLapackeEvd<std::complex<double>>;
} // namespace pairinteraction
