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

// Spherical entries as a dense matrix, entries that are not stored are zero
template <typename Scalar>
Eigen::MatrixX<Scalar>
get_constant_entries_as_matrix(const GreenTensorInterpolator<Scalar> &interpolator, int kappa1,
                               int kappa2) {
    Eigen::MatrixX<Scalar> matrix =
        Eigen::MatrixX<Scalar>::Zero(spherical::get_transformator<Scalar>(kappa1).rows(),
                                     spherical::get_transformator<Scalar>(kappa2).rows());
    for (const auto &[key, value] : get_constant_entries_as_map(interpolator, kappa1, kappa2)) {
        matrix(key.first, key.second) = value;
    }
    return matrix;
}

// Combinations of multipole orders for which the green tensors are implemented
constexpr std::array<std::pair<int, int>, 9> IMPLEMENTED_KAPPAS{
    {{0, 0}, {0, 1}, {1, 0}, {0, 2}, {2, 0}, {1, 1}, {1, 2}, {2, 1}, {2, 2}}};

// Polynomial in the cartesian components x, y, z, stored as a map from exponents to coefficients
using Polynomial = std::map<std::array<int, 3>, double>;

// Derivative of P(x) / r^(2p+1) with respect to x_a, returned as P'(x) / r^(2p+3) with
// P' = (d_a P) r^2 - (2p+1) x_a P
Polynomial differentiate(const Polynomial &poly, int p, int a) {
    Polynomial result;
    for (const auto &[exponents, coefficient] : poly) {
        if (exponents.at(a) > 0) {
            auto e = exponents;
            e.at(a) -= 1;
            for (int b = 0; b < 3; ++b) {
                auto e2 = e;
                e2.at(b) += 2;
                result[e2] += coefficient * exponents.at(a);
            }
        }
        auto e = exponents;
        e.at(a) += 1;
        result[e] -= (2 * p + 1) * coefficient;
    }
    return result;
}

// Reference for the cartesian green tensor of the multipole expansion,
// (-1)^kappa1 * d_{i_1} ... d_{i_{kappa1+kappa2}} (1/R), obtained by differentiating the Coulomb
// interaction symbolically. Rows (columns) enumerate the first kappa1 (last kappa2) indices with
// the first index being the least significant one.
Eigen::MatrixXd reference_cartesian_tensor(int kappa1, int kappa2,
                                           const Eigen::Vector3d &distance_vector) {
    const int order = kappa1 + kappa2;
    const double distance = distance_vector.norm();
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
                value += coefficient * std::pow(distance_vector.x(), exponents[0]) *
                    std::pow(distance_vector.y(), exponents[1]) *
                    std::pow(distance_vector.z(), exponents[2]);
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
}

DOCTEST_TEST_CASE("multipole green tensors for an arbitrary axis") {
    // The cartesian tensors are obtained by symbolic differentiation of 1/R and transformed to the
    // spherical basis with the same transformators as the implementation. This tests the cartesian
    // interaction tensors of GreenTensorInterpolator.cpp in isolation, errors in the transformators
    // cancel and are covered by spherical.test.cpp instead. The derivatives are symmetric under
    // permutations of the cartesian indices, so the order in which reference_cartesian_tensor()
    // decodes the row and column into indices does not matter.
    using complex_t = std::complex<double>;
    const std::vector<Eigen::Vector3d> distance_vectors{{0.3, -1.1, 0.7}};

    for (const auto &distance_vector : distance_vectors) {
        DOCTEST_INFO("distance_vector = (", distance_vector.x(), ", ", distance_vector.y(), ", ",
                     distance_vector.z(), ")");
        auto interpolator = GreenTensorInterpolator<complex_t>::from_multipole_expansion(
            {distance_vector.x(), distance_vector.y(), distance_vector.z()}, 5, 0, 0);

        for (const auto &kappas : IMPLEMENTED_KAPPAS) {
            // Structured bindings cannot be captured by the lambda inside DOCTEST_INFO before C++20
            const int kappa1 = kappas.first;
            const int kappa2 = kappas.second;
            DOCTEST_INFO("kappa1 = ", kappa1, ", kappa2 = ", kappa2);
            const Eigen::MatrixXcd actual =
                get_constant_entries_as_matrix(interpolator, kappa1, kappa2);

            const Eigen::MatrixXcd expected_cartesian =
                spherical::get_transformator<complex_t>(kappa1) *
                reference_cartesian_tensor(kappa1, kappa2, distance_vector).cast<complex_t>() *
                spherical::get_transformator<complex_t>(kappa2).adjoint();
            DOCTEST_REQUIRE(expected_cartesian.norm() > 0);
            DOCTEST_CHECK_MESSAGE((actual - expected_cartesian).norm() <
                                      1e-12 * expected_cartesian.norm(),
                                  "cartesian tensors");
        }
    }
}
} // namespace pairinteraction
