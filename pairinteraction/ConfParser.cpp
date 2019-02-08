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
