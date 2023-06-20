/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef SERIALIZATION_PATH_H
#define SERIALIZATION_PATH_H

#include "filesystem.hpp"
#include <cereal/types/string.hpp>

namespace cereal {

template <class Archive>
void serialize(Archive &ar, fs::path &p, unsigned int /* version */) {
    std::string s;

    if (Archive::is_saving::value) {
        s = p.string();
    }

    ar &cereal::make_nvp("path", s);

    if (Archive::is_loading::value) {
        p = fs::path(s);
    }
}

} // namespace cereal

#endif // SERIALIZATION_PATH_H
