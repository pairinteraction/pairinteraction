// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/Dense>
#include <array>
#include <complex>
#include <limits>
#include <stdexcept>

using namespace std::complex_literals;

namespace pairinteraction::spherical {

extern const Eigen::Matrix3<std::complex<double>> CARTESIAN_TO_SPHERICAL_KAPPA1;
extern const Eigen::Matrix<std::complex<double>, 6, 9> CARTESIAN_TO_SPHERICAL_KAPPA2;

template <typename Complex>
inline const Eigen::MatrixX<Complex> &get_transformator(int kappa) {
    static_assert(traits::NumTraits<Complex>::is_complex_v);
    if (kappa == 1) {
        static const auto mat = Eigen::MatrixX<Complex>(
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA1.template cast<Complex>());
        return mat;
    }
    if (kappa == 2) {
        static const auto mat = Eigen::MatrixX<Complex>(
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA2.template cast<Complex>());
        return mat;
    }
    throw std::invalid_argument("Invalid kappa value. Must be 1 or 2.");
}

template <typename Scalar>
inline std::array<Scalar, 3>
convert_to_spherical_basis(const std::array<typename traits::NumTraits<Scalar>::real_t, 3> &field) {
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    std::array<Scalar, 3> field_spherical{};
    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        Eigen::Map<Eigen::Vector3<Scalar>>(field_spherical.data(), field_spherical.size()) =
            CARTESIAN_TO_SPHERICAL_KAPPA1.template cast<Scalar>() *
            Eigen::Map<const Eigen::Vector3<real_t>>(field.data(), field.size());
    } else {
        if (std::abs(field[1]) > numerical_precision) {
            throw std::invalid_argument(
                "The field must not have a y-component if the scalar type is real.");
        }
        Eigen::Map<Eigen::Vector3<Scalar>>(field_spherical.data(), field_spherical.size()) =
            CARTESIAN_TO_SPHERICAL_KAPPA1.real().template cast<Scalar>() *
            Eigen::Map<const Eigen::Vector3<real_t>>(field.data(), field.size());
    }

    return field_spherical;
}
} // namespace pairinteraction::spherical
