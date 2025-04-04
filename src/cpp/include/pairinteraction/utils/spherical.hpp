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

template <typename Scalar>
inline const Eigen::MatrixX<Scalar> &get_transformator(int kappa) {
    if (kappa == 1) {
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            static const auto mat = Eigen::MatrixX<Scalar>(
                spherical::CARTESIAN_TO_SPHERICAL_KAPPA1.template cast<Scalar>());
            return mat;
        }
        static const auto mat = Eigen::MatrixX<Scalar>(
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA1.real().template cast<Scalar>());
        return mat;
    }
    if (kappa == 2) {
        if constexpr (traits::NumTraits<Scalar>::is_complex_v) {
            static const auto mat = Eigen::MatrixX<Scalar>(
                spherical::CARTESIAN_TO_SPHERICAL_KAPPA2.template cast<Scalar>());
            return mat;
        }
        static const auto mat = Eigen::MatrixX<Scalar>(
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA2.real().template cast<Scalar>());
        return mat;
    }
    throw std::invalid_argument("Invalid kappa value. Must be 1 or 2.");
}
} // namespace pairinteraction::spherical
