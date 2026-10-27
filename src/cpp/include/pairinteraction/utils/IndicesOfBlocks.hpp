// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <cstddef>

namespace pairinteraction {
struct IndicesOfBlock {
    IndicesOfBlock(size_t start, size_t end);
    size_t size() const;
    size_t start;
    size_t end;
};
} // namespace pairinteraction
