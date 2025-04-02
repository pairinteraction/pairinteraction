// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./Database.py.hpp"

#include "pairinteraction/basis/BasisAtom.hpp"
#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/enums/OperatorType.hpp"
#include "pairinteraction/ket/KetAtom.hpp"
#include "pairinteraction/operator/OperatorAtom.hpp"

#include <nanobind/eigen/sparse.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace pairinteraction;

static void declare_database(nb::module_ &m) {
    nb::class_<Database>(m, "Database")
        .def(nb::init<>())
        .def(nb::init<bool>(), "download_missing"_a)
        .def(nb::init<std::filesystem::path>(), "database_dir"_a)
        .def(nb::init<bool, bool, std::filesystem::path>(), "download_missing"_a, "use_cache"_a,
             "database_dir"_a);
}

void bind_database(nb::module_ &m) { declare_database(m); }
