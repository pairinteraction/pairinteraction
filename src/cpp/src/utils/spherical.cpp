// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/spherical.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <complex>

using namespace std::complex_literals;

namespace pairinteraction::spherical {

// The following matrices transform the ordered cartesian components of multipole operators
// (columns: x_a for kappa == 1, x_a*x_b for kappa == 2, x_a*x_b*x_c for kappa == 3, in row-major
// order) into the spherical multipole operators p_{kappa,q} = r^kappa * sqrt(4*pi/(2*kappa+1)) *
// Y_{kappa,q}. The rows labeled p_{kappa,q} contain the expansion of these operators in the
// cartesian monomials, divided by (2*kappa-1)!! = 1, 3, 15 for kappa = 1, 2, 3; the additional
// trace row p_{0,0} of the kappa == 2 matrix contains r^2/6. References:
// - https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
//   [explicit expressions for the spherical harmonics in cartesian coordinates]
// - D. A. Varshalovich, A. N. Moskalev, V. K. Khersonskii, Quantum Theory of Angular Momentum
//   (World Scientific, 1988), Chap. 5, https://doi.org/10.1142/0270
//   [tables of spherical and solid harmonics and their relation to cartesian tensors]
// - S. Weber et al., J. Phys. B 50, 133001 (2017), https://doi.org/10.1088/1361-6455/aa743a
//   [Eq. (8) defines the spherical multipole operators p_{kappa,q} used here]

const double SQRT_2 = std::sqrt(2.);
const double SQRT_2_3 = std::sqrt(2. / 3.);
const double SQRT_8_3 = std::sqrt(8. / 3.);
const double SQRT_3 = std::sqrt(3.);
const double SQRT_5 = std::sqrt(5.);
const double SQRT_30 = std::sqrt(30.);

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

// In contrast to the kappa == 2 case, no additional trace rows are included because the
// free-space green tensor is traceless in every pair of octupole indices; trace components of
// user-supplied cartesian tensors are projected out.
// clang-format off
// NOLINTNEXTLINE(cert-err58-cpp)
const Eigen::Matrix<std::complex<double>, 7, 27> CARTESIAN_TO_SPHERICAL_KAPPA3 =
(Eigen::Matrix<std::complex<double>, 7, 27>() <<
//  xxx           xxy           xxz           xyx           xyy           xyz           xzx           xzy           xzz           yxx           yxy           yxz           yyx           yyy           yyz           yzx           yzy           yzz           zxx           zxy           zxz           zyx           zyy           zyz           zzx           zzy           zzz
    3*SQRT_5,     -3i*SQRT_5,   0,            -3i*SQRT_5,   -3*SQRT_5,    0,            0,            0,            0,            -3i*SQRT_5,   -3*SQRT_5,    0,            -3*SQRT_5,    3i*SQRT_5,    0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            // p_{3,-3}
    0,            0,            SQRT_30,      0,            0,            -1i*SQRT_30,  SQRT_30,      -1i*SQRT_30,  0,            0,            0,            -1i*SQRT_30,  0,            0,            -SQRT_30,     -1i*SQRT_30,  -SQRT_30,     0,            SQRT_30,      -1i*SQRT_30,  0,            -1i*SQRT_30,  -SQRT_30,     0,            0,            0,            0,            // p_{3,-2}
    -3*SQRT_3,    1i*SQRT_3,    0,            1i*SQRT_3,    -SQRT_3,      0,            0,            0,            4*SQRT_3,     1i*SQRT_3,    -SQRT_3,      0,            -SQRT_3,      3i*SQRT_3,    0,            0,            0,            -4i*SQRT_3,   0,            0,            4*SQRT_3,     0,            0,            -4i*SQRT_3,   4*SQRT_3,     -4i*SQRT_3,   0,            // p_{3,-1}
    0,            0,            -6,           0,            0,            0,            -6,           0,            0,            0,            0,            0,            0,            0,            -6,           0,            -6,           0,            -6,           0,            0,            0,            -6,           0,            0,            0,            12,           // p_{3,0}
    3*SQRT_3,     1i*SQRT_3,    0,            1i*SQRT_3,    SQRT_3,       0,            0,            0,            -4*SQRT_3,    1i*SQRT_3,    SQRT_3,       0,            SQRT_3,       3i*SQRT_3,    0,            0,            0,            -4i*SQRT_3,   0,            0,            -4*SQRT_3,    0,            0,            -4i*SQRT_3,   -4*SQRT_3,    -4i*SQRT_3,   0,            // p_{3,1}
    0,            0,            SQRT_30,      0,            0,            1i*SQRT_30,   SQRT_30,      1i*SQRT_30,   0,            0,            0,            1i*SQRT_30,   0,            0,            -SQRT_30,     1i*SQRT_30,   -SQRT_30,     0,            SQRT_30,      1i*SQRT_30,   0,            1i*SQRT_30,   -SQRT_30,     0,            0,            0,            0,            // p_{3,2}
    -3*SQRT_5,    -3i*SQRT_5,   0,            -3i*SQRT_5,   3*SQRT_5,     0,            0,            0,            0,            -3i*SQRT_5,   3*SQRT_5,     0,            3*SQRT_5,     3i*SQRT_5,    0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0             // p_{3,3}
).finished() * (1./180.);
// clang-format on

} // namespace pairinteraction::spherical
