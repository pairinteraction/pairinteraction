#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Sparse>
#include <complex>
#include <map>
#include <unsupported/Eigen/Splines>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class GreenTensor {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    class Entry final {
    public:
        Entry(int row, int col, Eigen::Spline<real_t, 1> real_spline,
              Eigen::Spline<real_t, 1> imag_spline);
        Scalar val(double omega) const;
        int row() const noexcept;
        int col() const noexcept;

    private:
        int row_;
        int col_;
        Eigen::Spline<real_t, 1> real_spline;
        Eigen::Spline<real_t, 1> imag_spline;
    };

    void
    set_entries(int kappa1, int kappa2,
                const std::vector<Eigen::SparseMatrix<Scalar>> &tensors_in_cartesian_coordinates,
                const std::vector<double> &omegas);

    const std::vector<Entry> &get_entries(int kappa1, int kappa2) const;

private:
    std::map<std::pair<int, int>, std::vector<Entry>> entries_map;
};

extern template class GreenTensor<double>;
extern template class GreenTensor<std::complex<double>>;
} // namespace pairinteraction
