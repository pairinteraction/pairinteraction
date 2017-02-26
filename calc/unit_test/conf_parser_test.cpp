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

#include "ConfParser.h"
#define BOOST_TEST_MODULE Configuration parser test
#include <boost/test/unit_test.hpp>

#include <iostream>


BOOST_AUTO_TEST_CASE( conf_parser_test )
{
  Configuration c;

  // Test subscript setter
  c["first"]  << 2;
  c["second"] << 1.2;
  c["third"]  << "Hello World!";
  c["fourth"] << c["first"];
  c["third"] >> c["fifth"];

  // Test subscript getter and str() conversion
  BOOST_CHECK_EQUAL( c["first"] .str() , "2" );
  BOOST_CHECK_EQUAL( c["second"].str() , "1.2" );
  BOOST_CHECK_EQUAL( c["third"] .str() , "Hello World!" );
  BOOST_CHECK_EQUAL( c["fourth"].str() , "2" );
  BOOST_CHECK_EQUAL( c["fifth"] .str() , "Hello World!" );

  // Test subscript getter and type deduction
  int i;
  double d;
  std::string s;
  c["first"]  >> i;
  c["second"] >> d;
  c["third"]  >> s;
  BOOST_CHECK_EQUAL( i , 2 );
  BOOST_CHECK_EQUAL( d , 1.2 );
  BOOST_CHECK_EQUAL( s , "Hello World!" );

  // Test const methods
  Configuration const cc(c);
  BOOST_CHECK_EQUAL( cc["first"] .str() , "2" );
  BOOST_CHECK_EQUAL( cc["second"].str() , "1.2" );
  BOOST_CHECK_EQUAL( cc["third"] .str() , "Hello World!" );
  BOOST_CHECK_EQUAL( cc["fourth"].str() , "2" );
  BOOST_CHECK_EQUAL( cc["fifth"] .str() , "Hello World!" );
  BOOST_CHECK_THROW( cc["nonexistent"] , std::out_of_range );

  // Test comparison
  BOOST_CHECK( c == cc );

  // finder
  auto it = c.find("first");
  BOOST_CHECK_EQUAL( it->first        , "first" );
  BOOST_CHECK_EQUAL( it->second.str() , "2"     );

  auto cit = cc.find("first");
  BOOST_CHECK_EQUAL( cit->first        , "first" );
  BOOST_CHECK_EQUAL( cit->second.str() , "2"     );
}
