// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/spherical.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <complex>

using namespace std::complex_literals;

namespace pairinteraction::spherical {

const double SQRT_2 = std::sqrt(2.);
const double SQRT_2_3 = std::sqrt(2. / 3.);
const double SQRT_8_3 = std::sqrt(8. / 3.);

// clang-format off
// NOLINTNEXTLINE(cert-err58-cpp)
const Eigen::Matrix3<std::complex<double>> CARTESIAN_TO_SPHERICAL_KAPPA1 =
(Eigen::Matrix3<std::complex<double>>() <<
//  x           y           z
    1,          -1i,        0,      // p_{1,-1}
    0,          0,          SQRT_2, // p_{1,0}
    -1,         -1i,        0       // p_{1,1}
).finished() * std::sqrt(1./2.);
// clang-format on

// clang-format off
// NOLINTNEXTLINE(cert-err58-cpp)
const Eigen::Matrix<std::complex<double>, 6, 9> CARTESIAN_TO_SPHERICAL_KAPPA2 =
(Eigen::Matrix<std::complex<double>, 6, 9>() <<
// xx           xy          xz          yx          yy          yz          zx          zy          zz
    1,          -1i,        0,          -1i,        -1,         0,          0,          0,          0,          // p_{2,-2}
    0,          0,          1,          0,          0,          -1i,        1,          -1i,        0,          // p_{2,-1}
    -SQRT_2_3,  0,          0,          0,          -SQRT_2_3,  0,          0,          0,          SQRT_8_3,   // p_{2,0}
    0,          0,          -1,         0,          0,          -1i,        -1,         -1i,        0,          // p_{2,1}
    1,          1i,         0,          1i,         -1,         0,          0,          0,          0,          // p_{2,2}
    SQRT_2_3,   0,          0,          0,          SQRT_2_3,   0,          0,          0,          SQRT_2_3    // p_{0,0}
).finished() * std::sqrt(1./24.);
// clang-format on

} // namespace pairinteraction::spherical
