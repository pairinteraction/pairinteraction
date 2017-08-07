/*
 * Copyright (c) 2017 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "SQLite.h"
#include "QuantumDefect.h"
#include "Wavefunction.h"
#define BOOST_TEST_MODULE Wavefunctions test
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <iostream>

template < int l >
struct Fixture
{
  QuantumDefect qd;
  Fixture() : qd( "Rb", 79, l, 1.5 ) {};
};

typedef boost::mpl::list<Fixture<1>, Fixture<2>> Fixtures;

BOOST_FIXTURE_TEST_CASE_TEMPLATE(model_potentials, T, Fixtures, T)
{
  // There could be better coverage
  BOOST_CHECK( std::isnan( model_potential::V(T::qd, 0) ) );
  BOOST_CHECK( std::isnan( model_potential::g(T::qd, 0) ) );
}


BOOST_FIXTURE_TEST_CASE_TEMPLATE(numerovs_method, T, Fixtures, T)
{
  Numerov N( T::qd );
  auto const& x = N.axis();
  auto const& y = N.integrate();

  // Check for correct number of integration steps
  BOOST_CHECK_EQUAL( x.size(), 12087 );
  BOOST_CHECK_EQUAL( y.size(), 12087 );

  // Check for correct upper bound and decay to zero
  BOOST_CHECK( x.back() <= std::sqrt( 2*T::qd.n*(T::qd.n+15) ) );
  BOOST_CHECK_SMALL( y.back(), 1e-6 );
}


BOOST_FIXTURE_TEST_CASE_TEMPLATE(coulomb_functions, T, Fixtures, T)
{
  Whittaker W( T::qd );
  auto const& x = W.axis();
  auto const& y = W.integrate();

  // Check for correct number of integration steps
  BOOST_CHECK_EQUAL( x.size(), 12087 );
  BOOST_CHECK_EQUAL( y.size(), 12087 );

  // Check for correct upper bound and decay to zero
  BOOST_CHECK( x.back() <= 2*T::qd.n*(T::qd.n+15) );
  BOOST_CHECK_SMALL( y.back(), 1e-6 );
}


BOOST_FIXTURE_TEST_CASE_TEMPLATE(method_comparison, T, Fixtures, T)
{
  Numerov N( T::qd );
  Whittaker W( T::qd );
  auto const& nx = W.axis();
  auto const& ny = W.integrate();
  auto const& wx = W.axis();
  auto const& wy = W.integrate();

  // Check whether both have the same number of points
  BOOST_CHECK_EQUAL( nx.size(), wx.size() );
  size_t n = nx.size();

  // Compare pointwise
  for (size_t i = 0; i < n; ++i)
    BOOST_CHECK_CLOSE( ny[i], wy[i], 1e-16 );
}

BOOST_FIXTURE_TEST_CASE_TEMPLATE(integration, T, Fixtures, T)
{
  BOOST_CHECK_CLOSE(
    IntegrateRadialElement<Numerov  >(T::qd, 1, T::qd),
    IntegrateRadialElement<Whittaker>(T::qd, 1, T::qd),
    1e-3 ); // corresponds to 0.1% deviation
}
