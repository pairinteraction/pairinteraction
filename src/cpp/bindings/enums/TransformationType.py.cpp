// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./TransformationType.py.hpp"

#include "pairinteraction/enums/TransformationType.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_transformation_type(nb::module_ &m) {
    nb::enum_<TransformationType>(m, "TransformationType", nb::is_arithmetic())
        .value("IDENTITY", TransformationType::IDENTITY)
        .value("SORT_BY_KET", TransformationType::SORT_BY_KET)
        .value("SORT_BY_QUANTUM_NUMBER_F", TransformationType::SORT_BY_QUANTUM_NUMBER_F)
        .value("SORT_BY_QUANTUM_NUMBER_M", TransformationType::SORT_BY_QUANTUM_NUMBER_M)
        .value("SORT_BY_PARITY", TransformationType::SORT_BY_PARITY)
        .value("SORT_BY_ENERGY", TransformationType::SORT_BY_ENERGY)
        .value("ROTATE", TransformationType::ROTATE)
        .value("ARBITRARY", TransformationType::ARBITRARY);
}
