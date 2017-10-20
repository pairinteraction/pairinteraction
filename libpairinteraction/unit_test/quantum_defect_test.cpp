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

#include "EmbeddedDatabase.h"
#include "QuantumDefect.h"
#include "SQLite.h"
#define BOOST_TEST_MODULE Quantum defect test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(qd_test)
{
    QuantumDefect qd("Rb", 45, 1, 0.5);

    // Check whether the input was stored correctly
    BOOST_CHECK_EQUAL(qd.species, "Rb");
    BOOST_CHECK_EQUAL(qd.n, 45);
    BOOST_CHECK_EQUAL(qd.l, 1);
    BOOST_CHECK_EQUAL(qd.j, 0.5);

    // Check whether values are correctly read from the db
    EmbeddedDatabase db{};
    sqlite::statement stmt(db);
    stmt.set("select ac,Z,a1,a2,a3,a4,rc from model_potential where ( (element "
             "= 'Rb') and (L = 1) );");
    // The database should be consistent
    BOOST_CHECK_NO_THROW(stmt.prepare());
    BOOST_CHECK_NO_THROW(stmt.step());

    // Check the retrieved values
    double ac = stmt.get<double>(0);
    int Z = stmt.get<double>(1);
    double a1 = stmt.get<double>(2);
    double a2 = stmt.get<double>(3);
    double a3 = stmt.get<double>(4);
    double a4 = stmt.get<double>(5);
    double rc = stmt.get<double>(6);

    BOOST_CHECK_EQUAL(qd.ac, ac);
    BOOST_CHECK_EQUAL(qd.Z, Z);
    BOOST_CHECK_EQUAL(qd.a1, a1);
    BOOST_CHECK_EQUAL(qd.a2, a2);
    BOOST_CHECK_EQUAL(qd.a3, a3);
    BOOST_CHECK_EQUAL(qd.a4, a4);
    BOOST_CHECK_EQUAL(qd.rc, rc);
}
