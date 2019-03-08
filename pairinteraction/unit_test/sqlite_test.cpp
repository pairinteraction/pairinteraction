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
#include <sstream>
#include <string>
#include <vector>
#define BOOST_TEST_MODULE SQLite interface test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(sqlite_thread_safety) // NOLINT
{
    // The sqlite3_threadsafe() function returns zero if and only if
    // SQLite was compiled with mutexing code omitted due to the
    // SQLITE_THREADSAFE compile-time option being set to 0.
    BOOST_CHECK(sqlite3_threadsafe());
}

BOOST_AUTO_TEST_CASE(sqlite_query_test) // NOLINT
{
    BOOST_CHECK_THROW(sqlite::handle db("no such database", SQLITE_OPEN_READWRITE), sqlite::error);
    // Open an in-memory database for tests
    sqlite::handle db(":memory:");

    sqlite::statement stmt(db);
    BOOST_CHECK_NO_THROW(stmt.set("This is not valid SQL"));
    BOOST_CHECK_THROW(stmt.step(), sqlite::error);
    BOOST_CHECK_THROW(stmt.prepare(), sqlite::error);
    BOOST_CHECK_THROW(stmt.exec("Neither is this"), sqlite::error);

    // Check string calling
    std::string string_query("create table test(text,integer,real);");
    BOOST_CHECK_NO_THROW(stmt.reset());
    BOOST_CHECK_NO_THROW(stmt.set(string_query));
    BOOST_CHECK_NO_THROW(stmt.prepare());
    BOOST_CHECK_EQUAL(stmt.step(), false);

    // Check stringstream calling
    std::stringstream ss_query;
    ss_query << "insert into test values(?1,?2,?3);";
    BOOST_CHECK_NO_THROW(stmt.reset());
    BOOST_CHECK_NO_THROW(stmt.set(ss_query));
    // Insert some stuff
    BOOST_CHECK_NO_THROW(stmt.prepare());
    BOOST_CHECK_NO_THROW(stmt.bind(1, "Hello World!"));
    BOOST_CHECK_NO_THROW(stmt.bind(2, 1729));
    BOOST_CHECK_NO_THROW(stmt.bind(3, 0.5));
    BOOST_CHECK_EQUAL(stmt.step(), false);
    // Reuse the query set above
    BOOST_CHECK_NO_THROW(stmt.reset());
    BOOST_CHECK_NO_THROW(stmt.prepare());
    BOOST_CHECK_NO_THROW(stmt.bind(1, "Goodbye Earth!"));
    BOOST_CHECK_NO_THROW(stmt.bind(2, 42));
    BOOST_CHECK_NO_THROW(stmt.bind(3, 1.125));
    BOOST_CHECK_EQUAL(stmt.step(), false);
    BOOST_CHECK_THROW(stmt.step(), sqlite::error);

    // Check result
    BOOST_CHECK_NO_THROW(stmt.reset());
    BOOST_CHECK_NO_THROW(stmt.set("select * from test;"));
    BOOST_CHECK_NO_THROW(stmt.prepare());
    BOOST_CHECK_EQUAL(stmt.step(), true);
    BOOST_CHECK_EQUAL(stmt.get<std::string>(0), "Hello World!");
    BOOST_CHECK_EQUAL(stmt.get<int>(1), 1729);
    BOOST_CHECK_EQUAL(stmt.get<double>(2), 0.5);
    BOOST_CHECK_EQUAL(stmt.step(), true);
    BOOST_CHECK_EQUAL(stmt.get<std::string>(0), "Goodbye Earth!");
    BOOST_CHECK_EQUAL(stmt.get<int>(1), 42);
    BOOST_CHECK_EQUAL(stmt.get<double>(2), 1.125);
    BOOST_CHECK_EQUAL(stmt.step(), false);

    // Check iteration
    BOOST_CHECK_NO_THROW(stmt.reset());
    BOOST_CHECK_NO_THROW(stmt.set("select * from test;"));
    BOOST_CHECK_NO_THROW(stmt.prepare());

    int count = 0;
    for (auto &&r : stmt) {
        BOOST_CHECK_NO_THROW(r.get<std::string>(0));
        BOOST_CHECK_NO_THROW(r.get<int>(1));
        BOOST_CHECK_NO_THROW(r.get<double>(2));
        ++count;
    }
    BOOST_CHECK_EQUAL(count, 2);

#ifndef NDEBUG
    BOOST_CHECK_THROW(*(stmt.end()), std::out_of_range);
#endif
}
