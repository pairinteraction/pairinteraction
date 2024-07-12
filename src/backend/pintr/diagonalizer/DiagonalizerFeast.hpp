#pragma once

#include "pintr/interfaces/DiagonalizerInterface.hpp"
#include "pintr/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pintr {
template <typename Scalar>
class DiagonalizerFeast : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    DiagonalizerFeast(int m0);
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              real_t min_eigenvalue, real_t max_eigenvalue,
                              int precision) const override;

private:
    int m0;
};

extern template class DiagonalizerFeast<float>;
extern template class DiagonalizerFeast<double>;
extern template class DiagonalizerFeast<std::complex<float>>;
extern template class DiagonalizerFeast<std::complex<double>>;
} // namespace pintr
