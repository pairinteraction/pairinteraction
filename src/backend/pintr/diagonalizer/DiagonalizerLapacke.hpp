#pragma once

#include "pintr/interfaces/DiagonalizerInterface.hpp"
#include "pintr/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pintr {
template <typename Scalar>
class DiagonalizerLapacke : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              Range<real_t> allowed_range_of_evals = Range<real_t>(),
                              int precision = 4) const override;
};

extern template class DiagonalizerLapacke<float>;
extern template class DiagonalizerLapacke<double>;
extern template class DiagonalizerLapacke<std::complex<float>>;
extern template class DiagonalizerLapacke<std::complex<double>>;
} // namespace pintr
