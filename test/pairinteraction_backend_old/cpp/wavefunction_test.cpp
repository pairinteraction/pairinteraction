/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "QuantumDefect.hpp"
#include "SQLite.hpp"
#include "Wavefunction.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <iostream>

template <int l>
struct Fixture {
    QuantumDefect qd;
    Fixture() : qd("Rb", 79, l, 1.5){};
};

TEST_CASE_TEMPLATE("model_potentials", T, Fixture<1>, Fixture<2>) // NOLINT
{
    T const fixture;
    auto const qd = fixture.qd;
    // There could be better coverage
    CHECK(std::isnan(model_potential::V(qd, 0)));
    CHECK(std::isnan(model_potential::g(qd, 0)));
}

TEST_CASE_TEMPLATE("numerovs_method", T, Fixture<1>, Fixture<2>) // NOLINT
{
    T const fixture;
    auto const qd = fixture.qd;
    Numerov N(qd);
    auto const &xy = N.integrate();

    // Check for correct number of integration steps
    CHECK(xy.rows() == 12087);

    // Check for correct upper bound and decay to zero
    CHECK(xy(xy.rows() - 1, 0) <= std::sqrt(2 * qd.n * (qd.n + 15)));
    CHECK(xy(xy.rows() - 1, 1) == doctest::Approx(0.0).epsilon(1e-6));
}

#ifdef WITH_GSL
TEST_CASE_TEMPLATE("coulomb_functions", T, Fixture<1>, Fixture<2>) // NOLINT
{
    T const fixture;
    auto const qd = fixture.qd;
    Whittaker W(qd);
    auto const &xy = W.integrate();

    // Check for correct number of integration steps
    CHECK(xy.rows() == 12087);

    // Check for correct upper bound and decay to zero
    CHECK(xy(xy.rows() - 1, 0) <= 2 * qd.n * (qd.n + 15));
    CHECK(xy(xy.rows() - 1, 1) == doctest::Approx(0.0).epsilon(1e-6));
}

TEST_CASE_TEMPLATE("method_comparison", T, Fixture<1>, Fixture<2>) // NOLINT
{
    T const fixture;
    auto const qd = fixture.qd;
    Numerov N(qd);
    Whittaker W(qd);
    auto const &nxy = N.integrate();
    auto const &wxy = W.integrate();

    // Check whether both have the same number of points
    CHECK(nxy.rows() == wxy.rows());
    size_t n = nxy.rows();

    // Compare pointwise
    for (size_t i = 0; i < n; ++i) {
        CHECK(std::sqrt(nxy(i, 0)) * nxy(i, 1) - wxy(i, 1) == doctest::Approx(0.0).epsilon(1e-2));
    }
}
TEST_CASE_TEMPLATE("integration", T, Fixture<1>, Fixture<2>) // NOLINT
{
    T const fixture;
    auto const qd = fixture.qd;
    double mu_n = IntegrateRadialElement<Numerov>(qd, 1, qd);
    double mu_w = IntegrateRadialElement<Whittaker>(qd, 1, qd);
    CHECK(mu_n == doctest::Approx(mu_w).scale(1e-3)); // corresponds to 0.1% deviation
}
#endif
