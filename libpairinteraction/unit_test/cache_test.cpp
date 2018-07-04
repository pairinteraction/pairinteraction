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

#include "Cache.h"
#define BOOST_TEST_MODULE Configuration parser test
#include <boost/test/unit_test.hpp>
#include <thread>
#include <vector>

BOOST_AUTO_TEST_CASE(cache_test)
{
    Cache<std::string, int> cache;

    BOOST_CHECK_NO_THROW(cache.save("Hello world!", 1));

    boost::optional<int> oe;
    BOOST_CHECK_NO_THROW(oe = cache.restore("Hello world!"));
    BOOST_CHECK_EQUAL(oe.get(), 1);

    BOOST_CHECK_NO_THROW(cache.clear());
}

BOOST_AUTO_TEST_CASE(smash_test)
{
    Cache<std::string, int> cache;

    BOOST_CHECK_NO_THROW(cache.save("Hello world!", 1));
    BOOST_CHECK_THROW(cache.save("Hello world!", 2), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(parallel_cache_test)
{
    Cache<int, std::string> cache;

    std::vector<std::thread> threads(10);

    for (std::size_t i = 0; i < 10; ++i) {
        threads[i] = std::thread([&cache, i]() {
            cache.save(i, "Hello from thread " + std::to_string(i));
        });
}

    for (std::size_t i = 0; i < 10; ++i) {
        threads[i].join();
}

    for (std::size_t i = 0; i < 10; ++i) {
        BOOST_CHECK_EQUAL(cache.restore(i).get(),
                          "Hello from thread " + std::to_string(i));
}
}
