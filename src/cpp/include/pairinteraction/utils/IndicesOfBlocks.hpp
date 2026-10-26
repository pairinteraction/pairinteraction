// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <cstddef>
#include <initializer_list>
#include <set>
#include <vector>

namespace pairinteraction {
struct IndicesOfBlock {
    IndicesOfBlock(size_t start, size_t end);
    size_t size() const;
    size_t start;
    size_t end;
};

class IndicesOfBlocksCreator {
public:
    IndicesOfBlocksCreator(std::initializer_list<size_t> boundaries);
    void add(size_t boundary);
    std::vector<IndicesOfBlock> create() const;
    size_t size() const;

private:
    std::set<size_t> boundaries;
};
} // namespace pairinteraction
