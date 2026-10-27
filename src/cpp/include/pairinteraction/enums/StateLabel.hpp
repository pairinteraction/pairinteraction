// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace pairinteraction {
// The quantities by which the states of a basis can be labeled, sorted and, except for the energy,
// grouped into blocks
enum class StateLabel : unsigned char {
    QUANTUM_NUMBER_F,
    QUANTUM_NUMBER_M,
    PARITY,
    PRODUCT_OF_PARITIES,
    PARITY_UNDER_INVERSION,
    PARITY_UNDER_PERMUTATION,
    ENERGY
};
} // namespace pairinteraction
