// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace pairinteraction {
// The quantum numbers by which the states of a basis can be labeled, sorted and grouped into
// blocks
enum class StateLabel : unsigned char {
    QUANTUM_NUMBER_F,
    QUANTUM_NUMBER_M,
    PARITY,
    PRODUCT_OF_PARITIES,
    PARITY_UNDER_INVERSION,
    PARITY_UNDER_PERMUTATION
};
} // namespace pairinteraction
