// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/Info.hpp"

#include "./Info.py.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_info(nb::module_ &m) {
    nb::class_<Info>(m, "Info")
        .def_ro_static("with_mkl", &Info::with_mkl)
        .def_ro_static("has_eigen", &Info::has_eigen)
        .def_ro_static("has_lapacke_evd", &Info::has_lapacke_evd)
        .def_ro_static("has_lapacke_evr", &Info::has_lapacke_evr)
        .def_ro_static("has_feast", &Info::has_feast);
}
