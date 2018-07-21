/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef SERIALIZATION_PATH_H
#define SERIALIZATION_PATH_H

#include <boost/filesystem.hpp>

namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, boost::filesystem::path &p, const unsigned int version) {
    (void)version;

    std::string s;

    if (Archive::is_saving::value) {
        s = p.string();
    }

    ar &s;

    if (Archive::is_loading::value) {
        p = boost::filesystem::path(s);
    }
}

} // namespace serialization
} // namespace boost

#endif // SERIALIZATION_PATH_H
