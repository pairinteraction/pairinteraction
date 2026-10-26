// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/IndicesOfBlocks.hpp"
#include "pairinteraction/utils/eigen_assertion.hpp"

#include <Eigen/Core>
#include <vector>

namespace pairinteraction {
enum class SorterType : unsigned char;

class SorterBuilderInterface {
public:
    virtual ~SorterBuilderInterface() = default;
    virtual Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>
    get_sorter(const std::vector<SorterType> &labels) const = 0;
    virtual std::vector<IndicesOfBlock>
    get_indices_of_blocks(const std::vector<SorterType> &labels) const = 0;
};
} // namespace pairinteraction
