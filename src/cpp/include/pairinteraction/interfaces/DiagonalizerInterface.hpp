#pragma once

#include "pairinteraction/enums/FPP.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
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

    DiagonalizerInterface(FPP fpp);
    virtual ~DiagonalizerInterface() = default;
    virtual EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                      int precision) const = 0;
    virtual EigenSystemH<Scalar> eigh(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                      real_t min_eigenvalue, real_t max_eigenvalue,
                                      int precision) const;

protected:
    FPP fpp;
    template <typename ScalarLim>
    Eigen::MatrixX<ScalarLim> subtract_mean(const Eigen::MatrixX<Scalar> &matrix, real_t &shift,
                                            int precision) const;
    template <typename RealLim>
    Eigen::VectorX<real_t> add_mean(const Eigen::VectorX<RealLim> &eigenvalues, real_t shift) const;
};

extern template class DiagonalizerInterface<double>;
extern template class DiagonalizerInterface<std::complex<double>>;

extern template Eigen::MatrixX<float>
DiagonalizerInterface<double>::subtract_mean(const Eigen::MatrixX<double> &matrix, double &shift,
                                             int precision) const;
extern template Eigen::MatrixX<std::complex<float>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, int precision) const;

extern template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<float> &shifted_eigenvalues,
                                        double shift) const;
extern template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<float> &shifted_eigenvalues, double shift) const;

extern template Eigen::MatrixX<double>
DiagonalizerInterface<double>::subtract_mean(const Eigen::MatrixX<double> &matrix, double &shift,
                                             int precision) const;
extern template Eigen::MatrixX<std::complex<double>>
DiagonalizerInterface<std::complex<double>>::subtract_mean(
    const Eigen::MatrixX<std::complex<double>> &matrix, double &shift, int precision) const;

extern template Eigen::VectorX<double>
DiagonalizerInterface<double>::add_mean(const Eigen::VectorX<double> &shifted_eigenvalues,
                                        double shift) const;
extern template Eigen::VectorX<double> DiagonalizerInterface<std::complex<double>>::add_mean(
    const Eigen::VectorX<double> &shifted_eigenvalues, double shift) const;
} // namespace pairinteraction
