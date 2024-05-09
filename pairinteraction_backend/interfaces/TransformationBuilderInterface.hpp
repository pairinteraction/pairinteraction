#pragma once

#include "enums/TransformationType.hpp"
#include "utils/traits.hpp"

#include <Eigen/SparseCore>
#include <array>
#include <complex>
#include <vector>

template <typename Scalar>
struct Transformation {
    Eigen::SparseMatrix<Scalar, Eigen::RowMajor> matrix;
    std::vector<TransformationType> transformation_type;
};

struct Sorting {
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> matrix;
    std::vector<TransformationType> transformation_type;
};

struct Blocks {
    std::vector<int> start;
    std::vector<TransformationType> transformation_type;
};

template <typename Scalar>
class TransformationBuilderInterface {
public:
    using real_t = typename traits::NumTraits<Scalar>::real_t;

    virtual size_t get_number_of_states() const = 0;
    virtual size_t get_number_of_kets() const = 0;
    virtual const Transformation<Scalar> &get_transformation() const = 0;
    virtual Transformation<Scalar> get_rotator(real_t alpha, real_t beta, real_t gamma) const = 0;
    virtual Sorting get_sorter(TransformationType label) const = 0;
    virtual Blocks get_blocks(TransformationType label) const = 0;

    Transformation<Scalar> get_rotator(std::array<real_t, 3> to_z_axis,
                                       std::array<real_t, 3> to_y_axis) const;
};

extern template class TransformationBuilderInterface<float>;
extern template class TransformationBuilderInterface<double>;
extern template class TransformationBuilderInterface<std::complex<float>>;
extern template class TransformationBuilderInterface<std::complex<double>>;
