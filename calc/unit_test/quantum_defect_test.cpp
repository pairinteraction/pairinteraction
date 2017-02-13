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
#define BOOST_TEST_MODULE Quantum defect test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( qd_test )
{
  QuantumDefect qd( "Rb", 45, 1, 0.5 );

  // Check whether the input was stored correctly
  BOOST_CHECK_EQUAL( qd.species, "Rb" );
  BOOST_CHECK_EQUAL( qd.n      , 45   );
  BOOST_CHECK_EQUAL( qd.l      , 1    );
  BOOST_CHECK_EQUAL( qd.j      , 0.5  );

  // Check whether values are correctly read from the db
  sqlite::handle db("databases/quantum_defects.db");
  sqlite::result res = db.query( "select ac,Z,a1,a2,a3,a4,rc from model_potential "
                                 "where ( (element = 'Rb') and (L = 1) );" );
  // The database should be consistent
  BOOST_CHECK( res.size() > 0 );

  // Check the retrieved values
  int Z;
  double ac, a1, a2, a3, a4, rc;
  res.first() >> ac >> Z >> a1 >> a2 >> a3 >> a4 >> rc;

  BOOST_CHECK_EQUAL( qd.ac , ac );
  BOOST_CHECK_EQUAL( qd.Z  , Z  );
  BOOST_CHECK_EQUAL( qd.a1 , a1 );
  BOOST_CHECK_EQUAL( qd.a2 , a2 );
  BOOST_CHECK_EQUAL( qd.a3 , a3 );
  BOOST_CHECK_EQUAL( qd.a4 , a4 );
  BOOST_CHECK_EQUAL( qd.rc , rc );
}
