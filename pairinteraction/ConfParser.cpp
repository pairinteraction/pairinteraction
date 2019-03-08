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

#include "ConfParser.h"
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <string>

void Configuration::load_from_json(std::string const &filename) {
    using boost::property_tree::ptree;
    using boost::property_tree::json_parser::read_json;

    ptree pt;
    read_json(filename, pt);

    for (auto const &itr : pt) {
        params[itr.first] << itr.second.data();
    }
}

void Configuration::save_to_json(std::string const &filename) const {
    using boost::property_tree::ptree;
    using boost::property_tree::json_parser::write_json;

    ptree pt;

    for (auto const &p : *this) {
        pt.put(p.first, p.second.str());
    }

    write_json(filename, pt);
}
