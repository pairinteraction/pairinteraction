#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <complex>

namespace pairinteraction {
template <typename Scalar>
struct EigenSystemH {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> eigenvectors;
    Eigen::VectorX<real_t> eigenvalues;
};

template <typename Scalar>
class DiagonalizerInterface {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    virtual ~DiagonalizerInterface() = default;
    virtual EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                      int precision) const = 0;
    virtual EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                      real_t min_eigenvalue, real_t max_eigenvalue,
                                      int precision) const;
};

extern template class DiagonalizerInterface<double>;
extern template class DiagonalizerInterface<std::complex<double>>;
} // namespace pairinteraction
