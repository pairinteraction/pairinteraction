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
#include <sstream>
#include <string>
#include <vector>
#define BOOST_TEST_MODULE SQLite interface test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(sqlite_query_test)
{
    BOOST_CHECK_THROW(
        sqlite::handle db("no such database", SQLITE_OPEN_READWRITE),
        sqlite::error);
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
    BOOST_CHECK_THROW(stmt.set("whatever"), sqlite::error);
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
