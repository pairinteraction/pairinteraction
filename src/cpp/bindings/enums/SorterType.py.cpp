// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./SorterType.py.hpp"

#include "pairinteraction/enums/SorterType.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_sorter_type(nb::module_ &m) {
    nb::enum_<SorterType>(m, "SorterType")
        .value("QUANTUM_NUMBER_F", SorterType::QUANTUM_NUMBER_F)
        .value("QUANTUM_NUMBER_M", SorterType::QUANTUM_NUMBER_M)
        .value("PARITY", SorterType::PARITY)
        .value("PARITY_UNDER_INVERSION", SorterType::PARITY_UNDER_INVERSION)
        .value("PARITY_UNDER_PERMUTATION", SorterType::PARITY_UNDER_PERMUTATION)
        .value("ENERGY", SorterType::ENERGY);
}
