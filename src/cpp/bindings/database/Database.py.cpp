// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Database.py.hpp"

#include "pairinteraction/database/Database.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

static void declare_database(nb::module_ &m) {
    nb::class_<Database>(m, "Database")
        .def(nb::init<>())
        .def(nb::init<bool>(), "download_missing"_a)
        .def(nb::init<std::filesystem::path>(), "database_dir"_a)
        .def(nb::init<bool, bool, std::filesystem::path>(), "download_missing"_a, "use_cache"_a,
             "database_dir"_a)
        .def("get_download_missing", &Database::get_download_missing)
        .def("get_use_cache", &Database::get_use_cache)
        .def("get_database_dir", &Database::get_database_dir);
}

void bind_database(nb::module_ &m) { declare_database(m); }
