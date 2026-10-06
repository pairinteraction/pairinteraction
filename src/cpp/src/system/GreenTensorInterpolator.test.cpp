// SPDX-FileCopyrightText: 2025 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/GreenTensorInterpolator.hpp"

#include "pairinteraction/utils/spherical.hpp"

#include <Eigen/Dense>
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
} // namespace pairinteraction
