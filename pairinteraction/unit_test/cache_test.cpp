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

#include "Cache.h"
#define BOOST_TEST_MODULE Configuration parser test
#include <boost/test/unit_test.hpp>

#include <optional>
#include <thread>
#include <vector>

BOOST_AUTO_TEST_CASE(cache_test) // NOLINT
{
    Cache<std::string, int> cache;

    BOOST_CHECK_NO_THROW(cache.save("Hello world!", 1));

    std::optional<int> oe;
    BOOST_CHECK_NO_THROW(oe = cache.restore("Hello world!"));
    BOOST_CHECK_EQUAL(oe.value(), 1);

    BOOST_CHECK_NO_THROW(cache.clear());
}

BOOST_AUTO_TEST_CASE(smash_test) // NOLINT
{
    Cache<std::string, int> cache;

    BOOST_CHECK_NO_THROW(cache.save("Hello world!", 1));
    BOOST_CHECK_THROW(cache.save("Hello world!", 2), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(parallel_cache_test) // NOLINT
{
    Cache<int, std::string> cache;

    std::vector<std::thread> threads(10);

    for (std::size_t i = 0; i < 10; ++i) {
        threads[i] =
            std::thread([&cache, i]() { cache.save(i, "Hello from thread " + std::to_string(i)); });
    }

    for (std::size_t i = 0; i < 10; ++i) {
        threads[i].join();
    }

    for (std::size_t i = 0; i < 10; ++i) {
        BOOST_CHECK_EQUAL(cache.restore(i).value(), "Hello from thread " + std::to_string(i));
    }
}
