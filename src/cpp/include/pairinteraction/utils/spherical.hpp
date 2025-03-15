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
#include <unsupported/Eigen/KroneckerProduct>

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

template <typename Scalar, std::size_t Order>
inline std::array<Scalar, 3 * Order> convert_to_spherical_basis(
    const std::array<typename traits::NumTraits<Scalar>::real_t, 3> &vector) {
    static_assert(Order == 1 || Order == 2, "The order must be 1 or 2.");
    using real_t = typename traits::NumTraits<Scalar>::real_t;
    constexpr real_t numerical_precision = 100 * std::numeric_limits<real_t>::epsilon();

    Eigen::MatrixX<std::complex<double>> CARTESIAN_TO_SPHERICAL_KAPPA;
    if constexpr (Order == 1) {
        CARTESIAN_TO_SPHERICAL_KAPPA = CARTESIAN_TO_SPHERICAL_KAPPA1;
    } else {
        CARTESIAN_TO_SPHERICAL_KAPPA = CARTESIAN_TO_SPHERICAL_KAPPA2;
    }

    Eigen::MatrixX<Scalar> converter;
    if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
        converter = CARTESIAN_TO_SPHERICAL_KAPPA.template cast<Scalar>();
    } else {
        converter = CARTESIAN_TO_SPHERICAL_KAPPA.real().template cast<Scalar>();
        if (std::abs(vector[1]) > numerical_precision) {
            throw std::invalid_argument(
                "The vector must not have a y-component if the scalar type is real.");
        }
    }

    Eigen::Map<const Eigen::Vector3<real_t>> vector_map(vector.data(), vector.size());

    // Calculate sqrt(4pi/3) * r * Y_1,q with q = -1, 0, 1 for the first three elements
    std::array<Scalar, 3 * Order> vector_spherical{};
    if constexpr (Order == 1) {
        Eigen::Map<Eigen::Vector3<Scalar>>(vector_spherical.data(), vector_spherical.size()) =
            converter * vector_map;
        return vector_spherical;
    }

    // Calculate sqrt(4pi/5) * r^2 * Y_2,q / 3 with q = -2, -1, 0, 1, 2 for the first five elements
    // and (x^2 + y^2 + z^2) / 6 as the last element
    Eigen::Map<Eigen::Vector<Scalar, 6>>(vector_spherical.data(), vector_spherical.size()) =
        converter * Eigen::KroneckerProduct(vector_map, vector_map);
    return vector_spherical;
}
} // namespace pairinteraction::spherical
