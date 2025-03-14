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

    Eigen::Matrix3<Scalar> converter;
    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        converter = CARTESIAN_TO_SPHERICAL_KAPPA1;
    } else {
        converter = CARTESIAN_TO_SPHERICAL_KAPPA1.real();
        if (std::abs(field[1]) > numerical_precision) {
            throw std::invalid_argument(
                "The vector must not have a y-component if the scalar type is real.");
        }
    }

    std::array<Scalar, 3> field_spherical{};
    Eigen::Map<Eigen::Vector3<Scalar>>(field_spherical.data(), field_spherical.size()) =
        converter.template cast<Scalar>() *
        Eigen::Map<const Eigen::Vector3<real_t>>(field.data(), field.size());
    return field_spherical;
}

template <typename Scalar, int Order>
inline std::array<Scalar, 2 * Order + 1> get_multipole_expansion_factors(
    const std::array<typename traits::NumTraits<Scalar>::real_t, 3> &vector) {
    static_assert(Order == 1 || Order == 2, "The order must be 1 or 2.");
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    Scalar ii = 0;
    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        ii = Scalar(0, 1);
    } else {
        if (std::abs(vector[1]) > numerical_precision) {
            throw std::invalid_argument(
                "The vector must not have a y-component if the scalar type is real.");
        }
    }

    real_t norm = Eigen::Map<const Eigen::Vector3<real_t>>(vector.data(), vector.size()).norm();
    real_t scaling = 1 / std::pow(norm, Order + 1);
    real_t x = vector[0] / norm;
    real_t y = vector[1] / norm;
    real_t z = vector[2] / norm;

    // Calculate sqrt(4pi/3) * (-1)^q * Y_1,-q
    std::array<Scalar, 2 * Order + 1> factors{};
    if constexpr (Order == 1) {
        const real_t INV_SQRT_2 = 1 / std::sqrt(real_t(2));
        factors[0] = scaling * INV_SQRT_2 * (ii * y + x); // q = -1
        factors[1] = scaling * z;                         // q = 0
        factors[2] = scaling * INV_SQRT_2 * (ii * y - x); // q = 1
        return factors;
    }

    // Calculate sqrt(4pi/5) * (-1)^q * Y_2,-q
    const real_t SQRT_3_OVER_8 = std::sqrt(real_t(3) / real_t(8));
    const real_t SQRT_3_OVER_2 = std::sqrt(real_t(3) / real_t(2));
    factors[0] = scaling * SQRT_3_OVER_8 * ((ii * y + x) * (ii * y + x)); // q = -2
    factors[1] = scaling * SQRT_3_OVER_2 * z * (ii * y + x);              // q = -1
    factors[2] = scaling * 0.5 * (3 * z * z - 1);                         // q =  0
    // NOLINTBEGIN(cppcoreguidelines-pro-bounds-constant-array-index)
    factors[3] = scaling * SQRT_3_OVER_2 * z * (ii * y - x);              // q =  1
    factors[4] = scaling * SQRT_3_OVER_8 * ((ii * y - x) * (ii * y - x)); // q =  2
    // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index)
    return factors;
}
} // namespace pairinteraction::spherical
