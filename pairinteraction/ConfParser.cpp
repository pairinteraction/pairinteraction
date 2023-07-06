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

#include <cereal/external/rapidjson/document.h>
#include <cereal/external/rapidjson/istreamwrapper.h>
#include <cereal/external/rapidjson/ostreamwrapper.h>
#include <cereal/external/rapidjson/writer.h>

#include <fstream>
#include <string>

void Configuration::load_from_json(std::string const &filename) {
    std::ifstream ifs(filename);

    rapidjson::IStreamWrapper isw(ifs);
    rapidjson::Document document;
    document.ParseStream<>(isw);

    for (auto const &itr : document.GetObject()) {
        std::string name = itr.name.GetString();
        std::string value;
        switch (itr.value.GetType()) {
        case 0:
            // empty string
            break;
        case 1:
        case 2:
            value = itr.value.GetBool() ? "true" : "false";
            break;
        case 5:
            value = itr.value.GetString();
            break;
        case 6:
            value = std::to_string(itr.value.GetDouble());
            break;
        default:
            throw std::runtime_error("Can't convert value " + name);
        };
        params[name] << value;
    }
}

void Configuration::save_to_json(std::string const &filename) const {
    rapidjson::Document document;

    document.SetObject();
    for (auto const &p : *this) {
        rapidjson::Value name;
        name.SetString(p.first.c_str(), p.first.size(), document.GetAllocator());

        rapidjson::Value value;
        value.SetString(p.second.str().c_str(), p.second.str().size(), document.GetAllocator());

        document.AddMember(name, value, document.GetAllocator());
    }

    std::ofstream ofs(filename);

    rapidjson::OStreamWrapper osw(ofs);

    rapidjson::Writer<rapidjson::OStreamWrapper> writer(osw);
    document.Accept(writer);
}
