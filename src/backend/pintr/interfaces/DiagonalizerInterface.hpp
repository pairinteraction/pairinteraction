#pragma once

#include "pintr/utils/eigen_assertion.hpp"
#include "pintr/utils/eigen_compat.hpp"
#include "pintr/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <complex>
#include <semaphore>

namespace pintr {
// TODO replace this mocked class with the actual implementation
template <typename Real>
class Range {
public:
    bool is_finite() { return false; }
    Real min = -1000000;
    Real max = 1000000;
};

template <typename Scalar>
struct EigenSystemH {
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> evecs;
    Eigen::VectorX<real_t> evals;
};

template <typename Scalar>
class DiagonalizerInterface {
public:
    DiagonalizerInterface(int num_cpu_cores = -1);

    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    virtual ~DiagonalizerInterface() = default;
    virtual EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                      Range<real_t> allowed_range_of_evals = Range<real_t>(),
                                      int precision = 4) const = 0;

protected:
    mutable std::counting_semaphore<> smph_cpu_cores;
};

extern template class DiagonalizerInterface<float>;
extern template class DiagonalizerInterface<double>;
extern template class DiagonalizerInterface<std::complex<float>>;
extern template class DiagonalizerInterface<std::complex<double>>;
} // namespace pintr
