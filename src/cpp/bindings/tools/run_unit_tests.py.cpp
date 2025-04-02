// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./run_unit_tests.py.hpp"

#include "pairinteraction/tools/run_unit_tests.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/filesystem.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_run_unit_tests(nb::module_ &m) {
    m.def("run_unit_tests",
          [](bool download_missing, bool use_cache, std::filesystem::path database_dir) {
              return run_unit_tests(0, nullptr, download_missing, use_cache,
                                    std::move(database_dir));
          });
}
