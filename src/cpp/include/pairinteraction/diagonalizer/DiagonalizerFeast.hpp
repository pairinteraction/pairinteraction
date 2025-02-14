#pragma once

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/interfaces/DiagonalizerInterface.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <optional>

namespace pairinteraction {
template <typename Scalar>
class DiagonalizerFeast : public DiagonalizerInterface<Scalar> {
public:
    using typename DiagonalizerInterface<Scalar>::real_t;

    DiagonalizerFeast(int m0, FPP fpp = FPP::FLOAT64);
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              double atol) const override;
    EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                              std::optional<real_t> min_eigenvalue,
                              std::optional<real_t> max_eigenvalue, double atol) const override;

private:
    int m0;
    template <typename ScalarLim>
    EigenSystemH<Scalar> dispatch_eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                       real_t min_eigenvalue, real_t max_eigenvalue,
                                       double atol) const;
};

extern template class DiagonalizerFeast<double>;
extern template class DiagonalizerFeast<std::complex<double>>;
} // namespace pairinteraction
