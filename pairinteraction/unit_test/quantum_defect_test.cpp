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

#include "EmbeddedDatabase.h"
#include "QuantumDefect.h"
#include "SQLite.h"
#define BOOST_TEST_MODULE Quantum defect test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(qd_test) // NOLINT
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

BOOST_AUTO_TEST_CASE(qd_errors) // NOLINT
{
    BOOST_CHECK_THROW(QuantumDefect qd("nop", 0, 0, 0), std::exception);
    BOOST_CHECK_THROW(QuantumDefect qd("Rb", 0, 10000, 0), std::exception);

    try {
        QuantumDefect qd("nop", 0, 0, 0);
    } catch (std::exception const &e) {
        BOOST_CHECK(e.what());
    }

    try {
        QuantumDefect qd("Rb", 0, 10000, 0);
    } catch (std::exception const &e) {
        BOOST_CHECK(e.what());
    }
}
