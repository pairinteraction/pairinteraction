// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/IndicesOfBlocks.hpp"

namespace pairinteraction {

IndicesOfBlock::IndicesOfBlock(size_t start, size_t end) : start(start), end(end) {}

size_t IndicesOfBlock::size() const { return end - start; }

IndicesOfBlocksCreator::IndicesOfBlocksCreator(std::initializer_list<size_t> boundaries)
    : boundaries(boundaries) {}

void IndicesOfBlocksCreator::add(size_t boundary) { boundaries.insert(boundary); }

std::vector<IndicesOfBlock> IndicesOfBlocksCreator::create() const {
    std::vector<IndicesOfBlock> blocks;
    if (boundaries.empty()) {
        return blocks;
    }

    auto it = boundaries.begin();
    size_t start = *it++;

    while (it != boundaries.end()) {
        blocks.emplace_back(start, *it);
        start = *it++;
    }

    return blocks;
}

size_t IndicesOfBlocksCreator::size() const {
    return boundaries.empty() ? 0 : boundaries.size() - 1;
}
} // namespace pairinteraction
