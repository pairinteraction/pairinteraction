// SPDX-FileCopyrightText: 2025 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/system/GreenTensorInterpolator.hpp"

#include <cmath>
#include <doctest/doctest.h>
#include <map>
#include <utility>
#include <variant>

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
    auto interpolator = GreenTensorInterpolator<double>::from_multipole_expansion({0, 0, 1}, 5);

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
} // namespace pairinteraction
