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

#include "SQLite.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <sstream>
#include <string>
#include <vector>

TEST_CASE("sqlite_thread_safety") // NOLINT
{
    // The sqlite3_threadsafe() function returns zero if and only if
    // SQLite was compiled with mutexing code omitted due to the
    // SQLITE_THREADSAFE compile-time option being set to 0.
    CHECK(sqlite3_threadsafe());
}

TEST_CASE("sqlite_query_test") // NOLINT
{
    CHECK_THROWS_AS(sqlite::handle("no such database", SQLITE_OPEN_READWRITE), sqlite::error);
    // Open an in-memory database for tests
    sqlite::handle db(":memory:");

    sqlite::statement stmt(db);
    CHECK_NOTHROW(stmt.set("This is not valid SQL"));
    CHECK_THROWS_AS(stmt.step(), sqlite::error);
    CHECK_THROWS_AS(stmt.prepare(), sqlite::error);
    CHECK_THROWS_AS(stmt.exec("Neither is this"), sqlite::error);

    // Check string calling
    std::string string_query("create table test(text,integer,real);");
    CHECK_NOTHROW(stmt.reset());
    CHECK_NOTHROW(stmt.set(string_query));
    CHECK_NOTHROW(stmt.prepare());
    CHECK(stmt.step() == false);

    // Check stringstream calling
    std::stringstream ss_query;
    ss_query << "insert into test values(?1,?2,?3);";
    CHECK_NOTHROW(stmt.reset());
    CHECK_NOTHROW(stmt.set(ss_query));
    // Insert some stuff
    CHECK_NOTHROW(stmt.prepare());
    CHECK_NOTHROW(stmt.bind(1, "Hello World!"));
    CHECK_NOTHROW(stmt.bind(2, 1729));
    CHECK_NOTHROW(stmt.bind(3, 0.5));
    CHECK(stmt.step() == false);
    // Reuse the query set above
    CHECK_NOTHROW(stmt.reset());
    CHECK_NOTHROW(stmt.prepare());
    CHECK_NOTHROW(stmt.bind(1, "Goodbye Earth!"));
    CHECK_NOTHROW(stmt.bind(2, 42));
    CHECK_NOTHROW(stmt.bind(3, 1.125));
    CHECK(stmt.step() == false);
    CHECK_THROWS_AS(stmt.step(), sqlite::error);

    // Check result
    CHECK_NOTHROW(stmt.reset());
    CHECK_NOTHROW(stmt.set("select * from test;"));
    CHECK_NOTHROW(stmt.prepare());
    CHECK(stmt.step() == true);
    CHECK(stmt.get<std::string>(0) == "Hello World!");
    CHECK(stmt.get<int>(1) == 1729);
    CHECK(stmt.get<double>(2) == 0.5);
    CHECK(stmt.step() == true);
    CHECK(stmt.get<std::string>(0) == "Goodbye Earth!");
    CHECK(stmt.get<int>(1) == 42);
    CHECK(stmt.get<double>(2) == 1.125);
    CHECK(stmt.step() == false);

    // Check iteration
    CHECK_NOTHROW(stmt.reset());
    CHECK_NOTHROW(stmt.set("select * from test;"));
    CHECK_NOTHROW(stmt.prepare());

    int count = 0;
    for (auto &&r : stmt) {
        CHECK_NOTHROW(r.get<std::string>(0));
        CHECK_NOTHROW(r.get<int>(1));
        CHECK_NOTHROW(r.get<double>(2));
        ++count;
    }
    CHECK(count == 2);

#ifndef NDEBUG
    CHECK_THROWS_AS(*(stmt.end()), std::out_of_range);
#endif
}
