// SPDX-FileCopyrightText: 2025 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/GreenTensorInterpolator.hpp"

#include "pairinteraction/utils/spherical.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <doctest/doctest.h>
#include <map>
#include <utility>
#include <variant>
#include <vector>

namespace pairinteraction {

namespace {
template <typename Scalar>
std::map<std::pair<int, int>, Scalar>
get_constant_entries_as_map(const GreenTensorInterpolator<Scalar> &interpolator, int kappa1,
                            int kappa2) {
    std::map<std::pair<int, int>, Scalar> map;
    for (const auto &entry : interpolator.get_spherical_entries(kappa1, kappa2)) {
        const auto &constant_entry =
            std::get<typename GreenTensorInterpolator<Scalar>::ConstantEntry>(entry);
        map[{constant_entry.row(), constant_entry.col()}] = constant_entry.val();
    }
    return map;
}

double factorial(int n) { return std::tgamma(n + 1); }

double binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0;
    }
    return factorial(n) / (factorial(k) * factorial(n - k));
}

// Regular solid harmonic p_{n,m}(r) = sqrt(4 pi / (2n + 1)) r^n Y_{n,m}(r) with the Condon-Shortley
// phase, evaluated from its explicit polynomial form, e.g., p_{1,1}(r) = -(x + iy) / sqrt(2).
std::complex<double> solid_harmonic(int n, int m, const Eigen::Vector3d &r) {
    const std::complex<double> plus(-r.x() / 2, -r.y() / 2); // -(x + iy) / 2
    const std::complex<double> minus(r.x() / 2, -r.y() / 2); // (x - iy) / 2
    std::complex<double> sum = 0;
    for (int k = std::max(0, -m); 2 * k <= n - m; ++k) {
        sum += std::pow(plus, m + k) * std::pow(minus, k) * std::pow(r.z(), n - m - 2 * k) /
            (factorial(m + k) * factorial(k) * factorial(n - m - 2 * k));
    }
    return std::sqrt(factorial(n + m) * factorial(n - m)) * sum;
}
// Polynomial in the cartesian components x, y, z, stored as a map from exponents to coefficients
using Polynomial = std::map<std::array<int, 3>, double>;

// Derivative of P(x) / r^(2p+1) with respect to x_a, returned as P'(x) / r^(2p+3) with
// P' = (d_a P) r^2 - (2p+1) x_a P
Polynomial differentiate(const Polynomial &poly, int p, int a) {
    Polynomial result;
    for (const auto &[exponents, coefficient] : poly) {
        if (exponents[a] > 0) {
            auto e = exponents;
            e[a] -= 1;
            for (int b = 0; b < 3; ++b) {
                auto e2 = e;
                e2[b] += 2;
                result[e2] += coefficient * exponents[a];
            }
        }
        auto e = exponents;
        e[a] += 1;
        result[e] -= (2 * p + 1) * coefficient;
    }
    return result;
}

// Reference for the cartesian green tensor of the multipole expansion,
// (-1)^kappa1 * d_{i_1} ... d_{i_{kappa1+kappa2}} (1/R), obtained by differentiating the Coulomb
// interaction symbolically. Rows (columns) enumerate the first kappa1 (last kappa2) indices with
// the first index being the least significant one.
Eigen::MatrixXd reference_cartesian_tensor(int kappa1, int kappa2,
                                           const std::array<double, 3> &distance_vector) {
    const int order = kappa1 + kappa2;
    const double distance = std::sqrt(distance_vector[0] * distance_vector[0] +
                                      distance_vector[1] * distance_vector[1] +
                                      distance_vector[2] * distance_vector[2]);
    const int rows = static_cast<int>(std::pow(3, kappa1));
    const int cols = static_cast<int>(std::pow(3, kappa2));
    Eigen::MatrixXd tensor(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            // Collect the cartesian indices of the derivatives
            std::vector<int> indices;
            for (int r = row, i = 0; i < kappa1; ++i, r /= 3) {
                indices.push_back(r % 3);
            }
            for (int c = col, i = 0; i < kappa2; ++i, c /= 3) {
                indices.push_back(c % 3);
            }

            Polynomial poly{{{0, 0, 0}, 1}};
            for (int p = 0; p < order; ++p) {
                poly = differentiate(poly, p, indices[p]);
            }

            double value = 0;
            for (const auto &[exponents, coefficient] : poly) {
                value += coefficient * std::pow(distance_vector[0], exponents[0]) *
                    std::pow(distance_vector[1], exponents[1]) *
                    std::pow(distance_vector[2], exponents[2]);
            }
            tensor(row, col) =
                (kappa1 % 2 == 0 ? 1 : -1) * value / std::pow(distance, 2 * order + 1);
        }
    }
    return tensor;
}
} // namespace

DOCTEST_TEST_CASE("spherical entries of the multipole green tensors for a z-oriented axis") {
    // For a distance vector along z with unit distance, the spherical entries must reproduce the
    // known coefficients (-1)^(kappa2+q) * sqrt(binom(kappa1+kappa2, kappa1+q) *
    // binom(kappa1+kappa2, kappa2+q)) of the multipole expansion, see Eq. (7) of
    // S. Weber et al., J. Phys. B 50, 133001 (2017), https://doi.org/10.1088/1361-6455/aa743a.
    // The factor (-1)^q is not contained in Eq. (7) because the reference couples the operators
    // p_{kappa1,q} and p_{kappa2,-q} whereas here the first operator is conjugated so that both
    // operators carry the same q, using p^dagger_{kappa,q} = (-1)^q p_{kappa,-q}.
    // The spherical entries couple the operators p_{kappa1,q}^dagger and p_{kappa2,q} with
    // q = row - kappa1 and q = col - kappa2, respectively.
    auto interpolator =
        GreenTensorInterpolator<double>::from_multipole_expansion({0, 0, 1}, 5, 1, 1);

    DOCTEST_SUBCASE("dipole-dipole") {
        auto map = get_constant_entries_as_map(interpolator, 1, 1);
        DOCTEST_REQUIRE(map.size() == 3);
        DOCTEST_CHECK(map.at({0, 0}) == doctest::Approx(1));
        DOCTEST_CHECK(map.at({1, 1}) == doctest::Approx(-2));
        DOCTEST_CHECK(map.at({2, 2}) == doctest::Approx(1));
    }

    DOCTEST_SUBCASE("dipole-quadrupole") {
        auto map = get_constant_entries_as_map(interpolator, 1, 2);
        DOCTEST_REQUIRE(map.size() == 3);
        DOCTEST_CHECK(map.at({0, 1}) == doctest::Approx(-std::sqrt(3)));
        DOCTEST_CHECK(map.at({1, 2}) == doctest::Approx(3));
        DOCTEST_CHECK(map.at({2, 3}) == doctest::Approx(-std::sqrt(3)));
    }

    DOCTEST_SUBCASE("quadrupole-dipole") {
        auto map = get_constant_entries_as_map(interpolator, 2, 1);
        DOCTEST_REQUIRE(map.size() == 3);
        DOCTEST_CHECK(map.at({1, 0}) == doctest::Approx(std::sqrt(3)));
        DOCTEST_CHECK(map.at({2, 1}) == doctest::Approx(-3));
        DOCTEST_CHECK(map.at({3, 2}) == doctest::Approx(std::sqrt(3)));
    }

    DOCTEST_SUBCASE("quadrupole-quadrupole") {
        auto map = get_constant_entries_as_map(interpolator, 2, 2);
        DOCTEST_REQUIRE(map.size() == 5);
        DOCTEST_CHECK(map.at({0, 0}) == doctest::Approx(1));
        DOCTEST_CHECK(map.at({1, 1}) == doctest::Approx(-4));
        DOCTEST_CHECK(map.at({2, 2}) == doctest::Approx(6));
        DOCTEST_CHECK(map.at({3, 3}) == doctest::Approx(-4));
        DOCTEST_CHECK(map.at({4, 4}) == doctest::Approx(1));
    }

    DOCTEST_SUBCASE("dipole-octupole") {
        auto map = get_constant_entries_as_map(interpolator, 1, 3);
        DOCTEST_REQUIRE(map.size() == 3);
        DOCTEST_CHECK(map.at({0, 2}) == doctest::Approx(std::sqrt(6)));
        DOCTEST_CHECK(map.at({1, 3}) == doctest::Approx(-4));
        DOCTEST_CHECK(map.at({2, 4}) == doctest::Approx(std::sqrt(6)));
    }

    DOCTEST_SUBCASE("octupole-dipole") {
        auto map = get_constant_entries_as_map(interpolator, 3, 1);
        DOCTEST_REQUIRE(map.size() == 3);
        DOCTEST_CHECK(map.at({2, 0}) == doctest::Approx(std::sqrt(6)));
        DOCTEST_CHECK(map.at({3, 1}) == doctest::Approx(-4));
        DOCTEST_CHECK(map.at({4, 2}) == doctest::Approx(std::sqrt(6)));
    }
}

DOCTEST_TEST_CASE("cartesian multipole green tensors agree with derivatives of the Coulomb "
                  "interaction") {
    // Tests the cartesian interaction tensors of GreenTensorInterpolator.cpp in isolation: the
    // reference is obtained by symbolic differentiation of 1/R and transformed to the spherical
    // basis with the same transformators as the implementation, so errors in the transformators
    // cancel and are covered by spherical.test.cpp and by the full-chain test below instead.
    // The derivatives are symmetric under permutations of the cartesian indices, so the order in
    // which reference_cartesian_tensor() decodes the row and column into indices does not matter.
    using complex_t = std::complex<double>;
    const std::array<double, 3> distance_vector{0.3, -1.1, 0.7};
    auto interpolator =
        GreenTensorInterpolator<complex_t>::from_multipole_expansion(distance_vector, 5, 0, 0);

    const std::vector<std::pair<int, int>> kappas{{0, 0}, {0, 1}, {1, 0}, {0, 2}, {2, 0}, {1, 1},
                                                  {1, 2}, {2, 1}, {2, 2}, {1, 3}, {3, 1}};
    for (const auto &[kappa1, kappa2] : kappas) {
        DOCTEST_CAPTURE(kappa1);
        DOCTEST_CAPTURE(kappa2);

        Eigen::MatrixX<complex_t> expected = spherical::get_transformator<complex_t>(kappa1) *
            reference_cartesian_tensor(kappa1, kappa2, distance_vector).cast<complex_t>() *
            spherical::get_transformator<complex_t>(kappa2).adjoint();

        Eigen::MatrixX<complex_t> actual =
            Eigen::MatrixX<complex_t>::Zero(expected.rows(), expected.cols());
        for (const auto &[key, value] : get_constant_entries_as_map(interpolator, kappa1, kappa2)) {
            actual(key.first, key.second) = value;
        }

        DOCTEST_REQUIRE(expected.norm() > 0);
        DOCTEST_CHECK((actual - expected).norm() < 1e-12 * expected.norm());
    }
}

DOCTEST_TEST_CASE("full chain from cartesian tensors to spherical green tensor entries for an "
                  "arbitrary axis") {
    // The spherical entries for an arbitrary distance vector R = |R| u follow from the Legendre
    // generating function of 1/|R + r2 - r1| and the translation theorem of the solid harmonics,
    //   (-1)^(kappa2+q1) * sqrt(binom(n+M, kappa1-q1) * binom(n-M, kappa1+q1))
    //   * p_{n,M}(u)^* / |R|^(n+1)
    // with n = kappa1 + kappa2 and M = q2 - q1. For u = e_z, only M = 0 contributes and the
    // coefficients of the z-oriented test case are recovered. As the formula is independent of
    // the cartesian tensors and of the cartesian-to-spherical transformators, it tests the full
    // chain: every coefficient that enters the interaction, the normalization of the
    // transformators, and the complex phases of the M != 0 entries. Entries of the trace
    // row/column of the quadrupole transformator must vanish.
    const Eigen::Vector3d distance_vector(0.6, 0.8, 2.4); // |R| = 2.6, generic orientation
    const double distance = distance_vector.norm();
    const Eigen::Vector3d unitvec = distance_vector / distance;

    auto interpolator = GreenTensorInterpolator<std::complex<double>>::from_multipole_expansion(
        {distance_vector.x(), distance_vector.y(), distance_vector.z()}, 5, 0, 0);

    const std::vector<std::pair<int, int>> implemented_kappas{
        {0, 0}, {0, 1}, {1, 0}, {0, 2}, {2, 0}, {1, 1}, {1, 2}, {2, 1}, {2, 2}, {1, 3}, {3, 1}};

    for (const auto &[kappa1, kappa2] : implemented_kappas) {
        DOCTEST_INFO("kappa1 = ", kappa1, ", kappa2 = ", kappa2);
        const int n = kappa1 + kappa2;

        const auto rows = spherical::get_transformator<double>(kappa1).rows();
        const auto cols = spherical::get_transformator<double>(kappa2).rows();
        Eigen::MatrixXcd actual = Eigen::MatrixXcd::Zero(rows, cols);
        const auto map = get_constant_entries_as_map(interpolator, kappa1, kappa2);
        DOCTEST_REQUIRE(!map.empty());
        for (const auto &[key, val] : map) {
            actual(key.first, key.second) = val;
        }

        Eigen::MatrixXcd expected = Eigen::MatrixXcd::Zero(rows, cols);
        for (int q1 = -kappa1; q1 <= kappa1; ++q1) {
            for (int q2 = -kappa2; q2 <= kappa2; ++q2) {
                const int M = q2 - q1;
                if (std::abs(M) > n) {
                    continue;
                }
                const double sign = (kappa2 + q1) % 2 == 0 ? 1 : -1;
                expected(q1 + kappa1, q2 + kappa2) = sign *
                    std::sqrt(binomial(n + M, kappa1 - q1) * binomial(n - M, kappa1 + q1)) *
                    std::conj(solid_harmonic(n, M, unitvec)) / std::pow(distance, n + 1);
            }
        }

        DOCTEST_CHECK((actual - expected).norm() <= 1e-12 * expected.norm());
    }
}
} // namespace pairinteraction
