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
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <optional>
#include <thread>
#include <vector>

TEST_CASE("cache_test") // NOLINT
{
    Cache<std::string, int> cache;

    CHECK_NOTHROW(cache.save("Hello world!", 1));

    std::optional<int> oe;
    CHECK_NOTHROW(oe = cache.restore("Hello world!"));
    CHECK(oe.value() == 1);

    CHECK_NOTHROW(cache.clear());
}

TEST_CASE("smash_test") // NOLINT
{
    Cache<std::string, int> cache;

    CHECK_NOTHROW(cache.save("Hello world!", 1));
    CHECK_THROWS_AS(cache.save("Hello world!", 2), std::runtime_error);
}

TEST_CASE("parallel_cache_test") // NOLINT
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
        CHECK(cache.restore(i).value() == "Hello from thread " + std::to_string(i));
    }
}
