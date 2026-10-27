// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/IndicesOfBlocks.hpp"

namespace pairinteraction {

IndicesOfBlock::IndicesOfBlock(size_t start, size_t end) : start(start), end(end) {}

size_t IndicesOfBlock::size() const { return end - start; }

} // namespace pairinteraction
