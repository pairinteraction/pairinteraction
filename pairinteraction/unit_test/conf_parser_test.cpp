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

#include "ConfParser.hpp"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <iostream>

TEST_CASE("conf_parser_test") // NOLINT
{
    Configuration c;

    // Test subscript setter
    c["first"] << 2;
    c["second"] << 1.2;
    c["third"] << "Hello World!";
    c["fourth"] << c["first"];
    c["third"] >> c["fifth"];

    // Test subscript getter and str() conversion
    CHECK(c["first"].str() == "2");
    CHECK(c["second"].str() == "1.2");
    CHECK(c["third"].str() == "Hello World!");
    CHECK(c["fourth"].str() == "2");
    CHECK(c["fifth"].str() == "Hello World!");
    CHECK(c.count("first") == 1);

    // Test subscript getter and type deduction
    int i;
    double d;
    std::string s;
    c["first"] >> i;
    c["second"] >> d;
    c["third"] >> s;
    CHECK(i == 2);
    CHECK(d == 1.2);
    CHECK(s == "Hello World!");

    // Test const methods
    Configuration const cc(c);
    CHECK(cc["first"].str() == "2");
    CHECK(cc["second"].str() == "1.2");
    CHECK(cc["third"].str() == "Hello World!");
    CHECK(cc["fourth"].str() == "2");
    CHECK(cc["fifth"].str() == "Hello World!");
    CHECK_THROWS_AS(cc["nonexistent"], std::out_of_range);
    CHECK(cc.count("first") == 1);

    // Test comparison
    CHECK(c == cc);

    // finder
    auto it = c.find("first");
    CHECK(it->first == "first");
    CHECK(it->second.str() == "2");

    auto cit = cc.find("first");
    CHECK(cit->first == "first");
    CHECK(cit->second.str() == "2");

    // iterate
    for (auto it : c) {
        static_cast<void>(it);
    }
    for (auto cit : cc) {
        static_cast<void>(cit);
    }

    // Test fusion
    c += cc;
}
