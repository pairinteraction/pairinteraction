#pragma once

#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerLapacke : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    DiagonalizerLapacke();
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              int precision) const override;
};

extern template class DiagonalizerLapacke<double>;
extern template class DiagonalizerLapacke<std::complex<double>>;
} // namespace pairinteraction
