// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "./StateLabel.py.hpp"

#include "pairinteraction/enums/StateLabel.hpp"

#include <nanobind/nanobind.h>

namespace nb = nanobind;
using namespace pairinteraction;

void bind_state_label(nb::module_ &m) {
    nb::enum_<StateLabel>(m, "StateLabel")
        .value("QUANTUM_NUMBER_F", StateLabel::QUANTUM_NUMBER_F)
        .value("QUANTUM_NUMBER_M", StateLabel::QUANTUM_NUMBER_M)
        .value("PARITY", StateLabel::PARITY)
        .value("PRODUCT_OF_PARITIES", StateLabel::PRODUCT_OF_PARITIES)
        .value("PARITY_UNDER_INVERSION", StateLabel::PARITY_UNDER_INVERSION)
        .value("PARITY_UNDER_PERMUTATION", StateLabel::PARITY_UNDER_PERMUTATION);
}
