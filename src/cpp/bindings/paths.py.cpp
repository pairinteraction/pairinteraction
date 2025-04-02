// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./paths.py.hpp"

#include "pairinteraction/utils/paths.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_paths(nb::module_ &m) {
    m.def("get_cache_directory", &paths::get_cache_directory);
    m.def("get_config_directory", &paths::get_config_directory);
}
