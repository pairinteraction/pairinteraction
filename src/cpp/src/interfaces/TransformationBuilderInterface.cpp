// SPDX-FileCopyrightText: 2024 PairInteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/interfaces/TransformationBuilderInterface.hpp"

#include "pairinteraction/enums/TransformationType.hpp"

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

template <typename Scalar>
Transformation<Scalar>::Transformation(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix,
                                       std::vector<TransformationType> transformation_type)
    : matrix(std::move(matrix)), transformation_type(std::move(transformation_type)) {}

template <typename Scalar>
Transformation<Scalar>::Transformation(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix)
    : matrix(std::move(matrix)), transformation_type({TransformationType::ARBITRARY}) {}

Sorting::Sorting(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix,
                 std::vector<TransformationType> transformation_type)
    : matrix(std::move(matrix)), transformation_type(std::move(transformation_type)) {}

Sorting::Sorting(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix)
    : matrix(std::move(matrix)), transformation_type({TransformationType::ARBITRARY}) {}

// Explicit instantiations
// NOLINTBEGIN(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)
#define INSTANTIATE_TRANSFORMATION(SCALAR)                                                         \
    template struct Transformation<SCALAR>;                                                        \
    template class TransformationBuilderInterface<SCALAR>;
// NOLINTEND(bugprone-macro-parentheses, cppcoreguidelines-macro-usage)

INSTANTIATE_TRANSFORMATION(double)
INSTANTIATE_TRANSFORMATION(std::complex<double>)

#undef INSTANTIATE_TRANSFORMATION
} // namespace pairinteraction
