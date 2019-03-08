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

#include "QuantumDefect.h"
#include "SQLite.h"
#include "Wavefunction.h"
#define BOOST_TEST_MODULE Wavefunctions test
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>

template <int l>
struct Fixture {
    QuantumDefect qd;
    Fixture() : qd("Rb", 79, l, 1.5){};
};

typedef boost::mpl::list<Fixture<1>, Fixture<2>> Fixtures;

BOOST_FIXTURE_TEST_CASE_TEMPLATE(model_potentials, T, Fixtures, T) // NOLINT
{
    // There could be better coverage
    BOOST_CHECK(std::isnan(model_potential::V(T::qd, 0)));
    BOOST_CHECK(std::isnan(model_potential::g(T::qd, 0)));
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(numerovs_method, T, Fixtures, T) // NOLINT
{
    Numerov N(T::qd);
    auto const &xy = N.integrate();

    // Check for correct number of integration steps
    BOOST_CHECK_EQUAL(xy.rows(), 12087);

    // Check for correct upper bound and decay to zero
    BOOST_CHECK(xy(xy.rows() - 1, 0) <= std::sqrt(2 * T::qd.n * (T::qd.n + 15)));
    BOOST_CHECK_SMALL(xy(xy.rows() - 1, 1), 1e-6);
}

#ifdef WITH_GSL
BOOST_FIXTURE_TEST_CASE_TEMPLATE(coulomb_functions, T, Fixtures, T) // NOLINT
{
    Whittaker W(T::qd);
    auto const &xy = W.integrate();

    // Check for correct number of integration steps
    BOOST_CHECK_EQUAL(xy.rows(), 12087);

    // Check for correct upper bound and decay to zero
    BOOST_CHECK(xy(xy.rows() - 1, 0) <= 2 * T::qd.n * (T::qd.n + 15));
    BOOST_CHECK_SMALL(xy(xy.rows() - 1, 1), 1e-6);
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(method_comparison, T, Fixtures, T) // NOLINT
{
    Numerov N(T::qd);
    Whittaker W(T::qd);
    auto const &nxy = N.integrate();
    auto const &wxy = W.integrate();

    // Check whether both have the same number of points
    BOOST_CHECK_EQUAL(nxy.rows(), wxy.rows());
    size_t n = nxy.rows();

    // Compare pointwise
    for (size_t i = 0; i < n; ++i) {
        BOOST_CHECK_SMALL(std::sqrt(nxy(i, 0)) * nxy(i, 1) - wxy(i, 1), 1e-2);
    }
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(integration, T, Fixtures, T) // NOLINT
{
    double mu_n = IntegrateRadialElement<Numerov>(T::qd, 1, T::qd);
    double mu_w = IntegrateRadialElement<Whittaker>(T::qd, 1, T::qd);
    BOOST_CHECK_CLOSE(mu_n, mu_w, 1e-3); // corresponds to 0.1% deviation
}
#endif
