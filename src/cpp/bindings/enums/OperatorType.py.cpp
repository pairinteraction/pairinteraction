// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./OperatorType.py.hpp"

#include "pairinteraction/enums/OperatorType.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_operator_type(nb::module_ &m) {
    nb::enum_<OperatorType>(m, "OperatorType")
        .value("ZERO", OperatorType::ZERO)
        .value("ENERGY", OperatorType::ENERGY)
        .value("ELECTRIC_DIPOLE", OperatorType::ELECTRIC_DIPOLE)
        .value("ELECTRIC_QUADRUPOLE", OperatorType::ELECTRIC_QUADRUPOLE)
        .value("ELECTRIC_QUADRUPOLE_ZERO", OperatorType::ELECTRIC_QUADRUPOLE_ZERO)
        .value("ELECTRIC_OCTUPOLE", OperatorType::ELECTRIC_OCTUPOLE)
        .value("MAGNETIC_DIPOLE", OperatorType::MAGNETIC_DIPOLE)
        .value("ARBITRARY", OperatorType::ARBITRARY);
}
