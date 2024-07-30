#pragma once

#include "pintr/interfaces/DiagonalizerInterface.hpp"
#include "pintr/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pintr {
template <typename Scalar>
class DiagonalizerEigen : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              int precision) const override;
};

extern template class DiagonalizerEigen<float>;
extern template class DiagonalizerEigen<double>;
extern template class DiagonalizerEigen<std::complex<float>>;
extern template class DiagonalizerEigen<std::complex<double>>;
} // namespace pintr
