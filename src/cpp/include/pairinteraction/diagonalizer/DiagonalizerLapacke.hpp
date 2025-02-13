#pragma once

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerLapacke : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    DiagonalizerLapacke(FPP fpp = FPP::FLOAT64);
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              int precision) const override;

private:
    template <typename ScalarLim>
    EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                       int precision) const;
};

extern template class DiagonalizerLapacke<double>;
extern template class DiagonalizerLapacke<std::complex<double>>;
} // namespace pairinteraction
