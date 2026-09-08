// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/spherical.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <doctest/doctest.h>
#include <iterator>
#include <vector>

namespace pairinteraction {

namespace {
double factorial(int n) {
    double result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

// Associated Legendre polynomial P_l^m(x) without the Condon-Shortley phase, calculated with the
// standard recurrence relations
double associated_legendre(int l, int m, double x) {
    double p_curr = 1;
    for (int i = 1; i <= m; ++i) {
        p_curr *= (2 * i - 1) * std::sqrt(1 - x * x);
    }
    if (l == m) {
        return p_curr;
    }
    double p_prev = p_curr;
    p_curr *= (2 * m + 1) * x;
    for (int ll = m + 2; ll <= l; ++ll) {
        double p_next = ((2 * ll - 1) * x * p_curr - (ll + m - 1) * p_prev) / (ll - m);
        p_prev = p_curr;
        p_curr = p_next;
    }
    return p_curr;
}

// Spherical multipole operator p_{kappa,q} = r^kappa * sqrt(4*pi/(2*kappa+1)) * Y_{kappa,q} with
// the Condon-Shortley phase convention, evaluated at the position v
std::complex<double> evaluate_spherical_multipole(int kappa, int q, const Eigen::Vector3d &v) {
    int q_abs = std::abs(q);
    double r = v.norm();
    double phase = q_abs % 2 == 0 ? 1 : -1;
    std::complex<double> val = phase *
        std::sqrt(factorial(kappa - q_abs) / factorial(kappa + q_abs)) * std::pow(r, kappa) *
        associated_legendre(kappa, q_abs, v.z() / r) *
        std::exp(std::complex<double>(0, q_abs * std::atan2(v.y(), v.x())));
    if (q < 0) {
        val = phase * std::conj(val);
    }
    return val;
}

// Calculate the cartesian-to-spherical transformation matrix for the given kappa independently of
// the hardcoded matrices: fit the monomial coefficients of the spherical multipole operators
// divided by (2*kappa-1)!! (and of the additional trace row r^2/6 for kappa == 2) to their values
// at sample points and distribute them uniformly over the ordered cartesian components
Eigen::MatrixXcd calculate_transformator(int kappa) {
    // Enumerate the independent monomials x^i * y^j * z^k with i + j + k == kappa
    std::vector<std::array<int, 3>> monomials;
    for (int i = 0; i <= kappa; ++i) {
        for (int j = 0; j <= kappa - i; ++j) {
            monomials.push_back({i, j, kappa - i - j});
        }
    }
    auto num_monomials = static_cast<int>(monomials.size());

    double double_factorial = 1;
    for (int i = 3; i <= 2 * kappa - 1; i += 2) {
        double_factorial *= i;
    }

    int num_rows = kappa == 2 ? 2 * kappa + 2 : 2 * kappa + 1;
    int num_points = 10 * num_monomials;

    // Evaluate the monomials and the spherical multipole operators at deterministic sample points
    Eigen::MatrixXcd lhs(num_points, num_monomials);
    Eigen::MatrixXcd rhs(num_points, num_rows);
    for (int s = 0; s < num_points; ++s) {
        double radius = 0.8 + 0.05 * s;
        double theta = 0.5 + 0.1 * s;
        double phi = 0.3 + 0.7 * s;
        Eigen::Vector3d v{radius * std::sin(theta) * std::cos(phi),
                          radius * std::sin(theta) * std::sin(phi), radius * std::cos(theta)};
        for (int n = 0; n < num_monomials; ++n) {
            lhs(s, n) = std::pow(v.x(), monomials[n][0]) * std::pow(v.y(), monomials[n][1]) *
                std::pow(v.z(), monomials[n][2]);
        }
        for (int q = -kappa; q <= kappa; ++q) {
            rhs(s, q + kappa) = evaluate_spherical_multipole(kappa, q, v) / double_factorial;
        }
        if (kappa == 2) {
            rhs(s, num_rows - 1) = v.squaredNorm() / 6;
        }
    }
    Eigen::MatrixXcd coefficients = lhs.colPivHouseholderQr().solve(rhs);
    DOCTEST_REQUIRE((lhs * coefficients).isApprox(rhs, 1e-9));

    // Distribute the monomial coefficients uniformly over the ordered cartesian components
    Eigen::MatrixXcd transformator =
        Eigen::MatrixXcd::Zero(num_rows, static_cast<int>(std::pow(3, kappa)));
    for (int col = 0; col < transformator.cols(); ++col) {
        std::array<int, 3> exponents{};
        for (int i = 0, remainder = col; i < kappa; ++i, remainder /= 3) {
            ++exponents.at(remainder % 3);
        }
        auto n = static_cast<int>(std::distance(
            monomials.begin(), std::find(monomials.begin(), monomials.end(), exponents)));
        double multiplicity = factorial(kappa) /
            (factorial(exponents[0]) * factorial(exponents[1]) * factorial(exponents[2]));
        transformator.col(col) = coefficients.row(n).transpose() / multiplicity;
    }
    return transformator;
}
} // namespace

DOCTEST_TEST_CASE("compare the hardcoded transformation matrices to calculated ones") {
    DOCTEST_SUBCASE("kappa == 1") {
        DOCTEST_CHECK(calculate_transformator(1).isApprox(
            Eigen::MatrixXcd(spherical::CARTESIAN_TO_SPHERICAL_KAPPA1), 1e-9));
    }

    DOCTEST_SUBCASE("kappa == 2") {
        DOCTEST_CHECK(calculate_transformator(2).isApprox(
            Eigen::MatrixXcd(spherical::CARTESIAN_TO_SPHERICAL_KAPPA2), 1e-9));
    }

    DOCTEST_SUBCASE("kappa == 3") {
        DOCTEST_CHECK(calculate_transformator(3).isApprox(
            Eigen::MatrixXcd(spherical::CARTESIAN_TO_SPHERICAL_KAPPA3), 1e-9));
    }
}

DOCTEST_TEST_CASE("convert cartesian to spherical basis") {
    DOCTEST_SUBCASE("kappa == 1") {
        auto identity = spherical::CARTESIAN_TO_SPHERICAL_KAPPA1 *
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA1.adjoint();

        DOCTEST_CHECK(identity.isApprox(Eigen::Matrix3<double>::Identity(), 1e-9));
    }

    DOCTEST_SUBCASE("kappa == 2") {
        auto diagonal = spherical::CARTESIAN_TO_SPHERICAL_KAPPA2 *
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA2.adjoint();

        DOCTEST_CHECK(diagonal.isDiagonal(1e-9));
    }

    DOCTEST_SUBCASE("kappa == 3") {
        auto gram = spherical::CARTESIAN_TO_SPHERICAL_KAPPA3 *
            spherical::CARTESIAN_TO_SPHERICAL_KAPPA3.adjoint();

        Eigen::Matrix<std::complex<double>, 7, 7> expected =
            Eigen::Matrix<std::complex<double>, 7, 7>::Identity() / 90.;
        DOCTEST_CHECK(gram.isApprox(expected, 1e-9));
    }
}
} // namespace pairinteraction
