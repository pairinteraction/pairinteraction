// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

namespace pairinteraction {
enum class SorterType : unsigned char {
    QUANTUM_NUMBER_F,
    QUANTUM_NUMBER_M,
    PARITY,
    PARITY_UNDER_INVERSION,
    PARITY_UNDER_PERMUTATION,
    ENERGY
};
} // namespace pairinteraction
