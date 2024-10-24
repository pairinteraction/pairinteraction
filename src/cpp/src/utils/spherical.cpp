#include "pairinteraction/utils/spherical.hpp"

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <complex>

using namespace std::complex_literals;

namespace pairinteraction::spherical {

// clang-format off
const Eigen::Matrix3<std::complex<double>> CARTESIAN_TO_SPHERICAL1 =
(Eigen::Matrix3<std::complex<double>>() <<
//  x                       y                       z
    1 / std::sqrt(2.),     -1i/std::sqrt(2.),       0,
    0,                      0,                      1,
    -1 / std::sqrt(2.),    -1i/std::sqrt(2.),       0
).finished();
// clang-format on

// clang-format off
const Eigen::Matrix<std::complex<double>, 5, 9> CARTESIAN_TO_SPHERICAL2 =
(Eigen::Matrix<std::complex<double>, 5, 9>() <<
// xx                       xy                      xz                      yx                      yy                      yz                      zx                      zy                      zz
    std::sqrt(3./8.),       -1i*std::sqrt(3./8.),   0,                      -1i*std::sqrt(3./8.),   -std::sqrt(3./8.),      0,                      0,                      0,                      0,
    0,                      0,                      std::sqrt(3./8.),       0,                      0,                      -1i*std::sqrt(3./8.),   std::sqrt(3./8.),       -1i*std::sqrt(3./8.),   0,
    -std::sqrt(1./4.),      0,                      0,                      0,                      -std::sqrt(1./4.),      0,                      0,                      0,                      1,
    0,                      0,                      -std::sqrt(3./8.),      0,                      0,                      -1i*std::sqrt(3./8.),   -std::sqrt(3./8.),      -1i*std::sqrt(3./8.),   0,
    std::sqrt(3./8.),       1i*std::sqrt(3./8.),    0,                      1i*std::sqrt(3./8.),    -std::sqrt(3./8.),      0,                      0,                      0,                      0
).finished();
// clang-format on

} // namespace pairinteraction::spherical
