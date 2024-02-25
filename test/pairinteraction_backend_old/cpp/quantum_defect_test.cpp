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

#include "EmbeddedDatabase.hpp"
#include "QuantumDefect.hpp"
#include "SQLite.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

TEST_CASE("qd_test") // NOLINT
{
    QuantumDefect qd("Rb", 45, 1, 0.5);

    // Check whether the input was stored correctly
    CHECK(qd.species == "Rb");
    CHECK(qd.n == 45);
    CHECK(qd.l == 1);
    CHECK(qd.j == 0.5);

    // Check whether values are correctly read from the db
    EmbeddedDatabase db{};
    sqlite::statement stmt(db);
    stmt.set("select ac,Z,a1,a2,a3,a4,rc from model_potential where ( (element "
             "= 'Rb') and (L = 1) );");
    // The database should be consistent
    CHECK_NOTHROW(stmt.prepare());
    CHECK_NOTHROW(stmt.step());

    // Check the retrieved values
    double ac = stmt.get<double>(0);
    int Z = stmt.get<double>(1);
    double a1 = stmt.get<double>(2);
    double a2 = stmt.get<double>(3);
    double a3 = stmt.get<double>(4);
    double a4 = stmt.get<double>(5);
    double rc = stmt.get<double>(6);

    CHECK(qd.ac == ac);
    CHECK(qd.Z == Z);
    CHECK(qd.a1 == a1);
    CHECK(qd.a2 == a2);
    CHECK(qd.a3 == a3);
    CHECK(qd.a4 == a4);
    CHECK(qd.rc == rc);
}

TEST_CASE("qd_errors") // NOLINT
{
    CHECK_THROWS_AS(QuantumDefect("nop", 0, 0, 0), std::exception);
    CHECK_THROWS_AS(QuantumDefect("Rb", 0, 10000, 0), std::exception);

    try {
        QuantumDefect qd("nop", 0, 0, 0);
    } catch (std::exception const &e) {
        CHECK(e.what());
    }

    try {
        QuantumDefect qd("Rb", 0, 10000, 0);
    } catch (std::exception const &e) {
        CHECK(e.what());
    }
}
