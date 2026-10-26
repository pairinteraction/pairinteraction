// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./SorterBuilderInterface.py.hpp"

#include "pairinteraction/interfaces/SorterBuilderInterface.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

static void declare_sorter_builder_interface(nb::module_ &m) {
    nb::class_<SorterBuilderInterface> pyclass(m, "SorterBuilderInterface");
}

void bind_sorter_builder_interface(nb::module_ &m) { declare_sorter_builder_interface(m); }
