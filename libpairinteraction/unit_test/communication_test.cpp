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

#include "Communication.h"
#include <chrono>
#include <string>
#include <thread>
#define BOOST_TEST_MODULE ZeroMQ interface test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(exceptions_test) // NOLINT
{
    auto context = zmq::context();
    auto socket = context.socket(ZMQ_REQ);

    BOOST_CHECK_THROW(context.socket(-42), zmq::error);
    BOOST_CHECK_THROW(socket.bind("garbage"), zmq::error);
    BOOST_CHECK_THROW(socket.connect("garbage"), zmq::error);

    try {
        context.socket(-42);
    } catch (zmq::error const &e) {
        e.what();
    }
}

BOOST_AUTO_TEST_CASE(send_test) // NOLINT
{
    auto context = zmq::context();
    constexpr static char msg[] = "Into the void...";
    constexpr static size_t len = sizeof(msg);

    std::thread sender([&context]() {
        auto publisher = context.socket(ZMQ_PUB);
        publisher.bind("tcp://*:5555");

        std::this_thread::sleep_for(
            std::chrono::seconds{1}); // wait for the others to receive
        BOOST_CHECK_EQUAL(publisher.send(msg), len - 1);
    });

    std::thread receiver([&context]() {
        auto subscriber = context.socket(ZMQ_SUB);
        subscriber.connect("tcp://localhost:5555");
        BOOST_CHECK_EQUAL(zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE, nullptr, 0),
                          0);

        char buf[len];
        BOOST_CHECK_EQUAL(zmq_recv(subscriber, buf, len, 0), len);
        BOOST_CHECK_EQUAL(msg, buf);
    });

    sender.join();
    receiver.join();
}
