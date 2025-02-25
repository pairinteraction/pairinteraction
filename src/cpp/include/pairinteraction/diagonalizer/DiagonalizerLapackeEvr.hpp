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
                              double atol) const override;
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              std::optional<real_t> min_eigenvalue,
                              std::optional<real_t> max_eigenvalue, double atol) const override;

private:
    template <typename ScalarLim>
    EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                       std::optional<real_t> min_eigenvalue,
                                       std::optional<real_t> max_eigenvalue, double atol) const;
};

extern template class DiagonalizerLapackeEvr<double>;
extern template class DiagonalizerLapackeEvr<std::complex<double>>;
} // namespace pairinteraction
