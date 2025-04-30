// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <complex>
#include <map>
#include <unsupported/Eigen/Splines>
#include <variant>
#include <vector>

namespace pairinteraction {
template <typename Scalar>
class GreenTensor {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;
    using complex_t = std::complex<real_t>;

    class ConstantEntry {
        friend class GreenTensor;

    public:
        Scalar val() const;
        int row() const noexcept;
        int col() const noexcept;

    private:
        ConstantEntry(int row, int col, Scalar val);
        int row_;
        int col_;
        Scalar val_;
    };

    class OmegaDependentEntry {
        friend class GreenTensor;

    public:
        Scalar val(double omega) const;
        int row() const noexcept;
        int col() const noexcept;

    private:
        OmegaDependentEntry(int row, int col, Eigen::Spline<real_t, 1> real_spline,
                            Eigen::Spline<real_t, 1> imag_spline);
        int row_;
        int col_;
        Eigen::Spline<real_t, 1> real_spline;
        Eigen::Spline<real_t, 1> imag_spline;
    };

    using Entry = std::variant<ConstantEntry, OmegaDependentEntry>;

    void set_entries(int kappa1, int kappa2,
                     const Eigen::MatrixX<Scalar> &tensor_in_cartesian_coordinates);
    void set_entries(int kappa1, int kappa2,
                     const std::vector<Eigen::MatrixX<Scalar>> &tensors_in_cartesian_coordinates,
                     const std::vector<double> &omegas);
    const std::vector<Entry> &get_entries(int kappa1, int kappa2) const;

private:
    std::map<std::pair<int, int>, std::vector<Entry>> entries_map;
};

extern template class GreenTensor<double>;
extern template class GreenTensor<std::complex<double>>;
} // namespace pairinteraction
