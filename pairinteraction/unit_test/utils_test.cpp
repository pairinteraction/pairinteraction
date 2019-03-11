/*
 * Copyright (c) 2018 Sebastian Weber, Henri Menke. All rights reserved.
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

#include "utils.h"
#define BOOST_TEST_MODULE Configuration parser test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(randint_test) // NOLINT
{
    int a = utils::randint(0, 10);
    BOOST_CHECK(0 <= a && a <= 10);
}
