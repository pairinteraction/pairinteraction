// SPDX-FileCopyrightText: 2025 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./FloatType.py.hpp"

#include "pairinteraction/enums/FloatType.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_float_type(nb::module_ &m) {
    nb::enum_<FloatType>(m, "FloatType")
        .value("FLOAT32", FloatType::FLOAT32)
        .value("FLOAT64", FloatType::FLOAT64);
}
