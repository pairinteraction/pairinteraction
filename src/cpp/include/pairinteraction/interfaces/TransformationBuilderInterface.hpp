// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/traits.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <complex>
#include <initializer_list>
#include <set>
#include <vector>

namespace pairinteraction {
enum class TransformationType : unsigned char;

template <typename Scalar>
struct Transformation {
    Transformation() = default;
    Transformation(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix,
                   std::vector<TransformationType> transformation_type);
    Transformation(Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix);
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix;
    std::vector<TransformationType> transformation_type;
};

struct Sorting {
    Sorting() = default;
    Sorting(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix,
            std::vector<TransformationType> transformation_type);
    Sorting(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix);
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix;
    std::vector<TransformationType> transformation_type;
};

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

template <typename Scalar>
class TransformationBuilderInterface {
public:
    static_assert(traits::NumTraits<Scalar>::from_floating_point_v);

    using real_t = typename traits::NumTraits<Scalar>::real_t;

    virtual ~TransformationBuilderInterface() = default;
    virtual const Transformation<Scalar> &get_transformation() const = 0;
    virtual Transformation<Scalar> get_rotator(real_t alpha, real_t beta, real_t gamma) const = 0;
    virtual Sorting get_sorter(const std::vector<TransformationType> &labels) const = 0;
    virtual std::vector<IndicesOfBlock>
    get_indices_of_blocks(const std::vector<TransformationType> &labels) const = 0;

    Transformation<Scalar> get_rotator(const std::array<real_t, 3> &to_z_axis,
                                       const std::array<real_t, 3> &to_y_axis) const;
};

extern template struct Transformation<double>;
extern template struct Transformation<std::complex<double>>;

extern template class TransformationBuilderInterface<double>;
extern template class TransformationBuilderInterface<std::complex<double>>;
} // namespace pairinteraction
