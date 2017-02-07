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

#include "SQLite.h"
#define BOOST_TEST_MODULE SQLite interface test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( sqlite_query_test )
{
  // Open an in-memory database for tests
  sqlite::handle db(":memory:");
  
  BOOST_CHECK_THROW( db.query("This is not valid SQL"), sqlite::sqlite_error );
}
