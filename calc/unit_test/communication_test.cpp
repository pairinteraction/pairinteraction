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

#include "communication.h"
#include <string>
#define BOOST_TEST_MODULE ZeroMQ interface test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( exceptions_test )
{
  auto context = zmq::context();
  auto socket  = context.socket(ZMQ_REQ);
  
  BOOST_CHECK_THROW( context.socket ( -42 )       , zmq::error );
  BOOST_CHECK_THROW( socket .bind   ( "garbage" ) , zmq::error );
  BOOST_CHECK_THROW( socket .connect( "garbage" ) , zmq::error );
}

BOOST_AUTO_TEST_CASE( send_test )
{
  auto context    = zmq::context();
  auto publisher  = context.socket(ZMQ_PUB);
  publisher.bind("tcp://*:5555");

  std::string msg("Into the void...");
  
  BOOST_CHECK_EQUAL( publisher.send(msg.c_str()), msg.size() );
}
