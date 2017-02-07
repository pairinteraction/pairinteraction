/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#include <string>
#include <sstream>
#include <vector>
#include "SQLite.h"
#define BOOST_TEST_MODULE SQLite interface test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( sqlite_query_test )
{
  // Open an in-memory database for tests
  sqlite::handle db(":memory:");

  BOOST_CHECK_THROW( db.query("This is not valid SQL"),
                     sqlite::sqlite_error );
  BOOST_CHECK_THROW( db.exec("Neither is this"),
                     sqlite::sqlite_error );

  // Check string calling
  std::string string_query( "create table test(text);" );
  BOOST_CHECK_EQUAL( db.query( string_query ).size(), 0 );

  // Check stringstream calling
  std::stringstream ss_query;
  ss_query << "insert into test values(\"Hello World!\");";
  ss_query << "insert into test values(\"Goodbye Earth!\");";
  BOOST_CHECK_EQUAL( db.query( ss_query ).size(), 0 );

  // Check result
  sqlite::result res = db.query("select * from test;");
  BOOST_CHECK_EQUAL( res.size(), 2 );
  std::string entry = *res.begin();
  BOOST_CHECK_EQUAL( entry, "Hello World!" );

  // Check if result is iterable
  std::vector < std::string > values;
  values.reserve( res.size() );
  for ( auto const& row : res )
  {
    values.push_back( std::string(row) );
  }
  BOOST_CHECK_EQUAL( values[0], "Hello World!"   );
  BOOST_CHECK_EQUAL( values[1], "Goodbye Earth!" );
}
