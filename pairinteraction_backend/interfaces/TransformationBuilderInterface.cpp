#include "interfaces/TransformationBuilderInterface.hpp"
#include "enums/TransformationType.hpp"
#include "utils/euler.hpp"

template <typename Scalar>
Transformation<Scalar>::Transformation(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix,
                                       const std::vector<TransformationType> &transformation_type)
    : matrix(matrix), transformation_type(transformation_type) {}

template <typename Scalar>
Transformation<Scalar>::Transformation(const Eigen::SparseMatrix<Scalar, Eigen::RowMajor> &matrix)
    : matrix(matrix), transformation_type({TransformationType::ARBITRARY}) {}

Sorting::Sorting(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &matrix,
                 const std::vector<TransformationType> &transformation_type)
    : matrix(matrix), transformation_type(transformation_type) {}

Sorting::Sorting(const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> &matrix)
    : matrix(matrix), transformation_type({TransformationType::ARBITRARY}) {}

Blocks::Blocks(const std::vector<int> &start,
               const std::vector<TransformationType> &transformation_type)
    : start(start), transformation_type(transformation_type) {}

Blocks::Blocks(const std::vector<int> &start)
    : start(start), transformation_type({TransformationType::ARBITRARY}) {}

template <typename Scalar>
Transformation<Scalar>
TransformationBuilderInterface<Scalar>::get_rotator(std::array<real_t, 3> to_z_axis,
                                                    std::array<real_t, 3> to_y_axis) const {
    auto euler_zyz_angles = euler::get_euler_angles(to_z_axis, to_y_axis);
    return this->get_rotator(euler_zyz_angles[0], euler_zyz_angles[1], euler_zyz_angles[2]);
}

// Explicit instantiations
template class Transformation<float>;
template class Transformation<double>;
template class Transformation<std::complex<float>>;
template class Transformation<std::complex<double>>;

template class TransformationBuilderInterface<float>;
template class TransformationBuilderInterface<double>;
template class TransformationBuilderInterface<std::complex<float>>;
template class TransformationBuilderInterface<std::complex<double>>;
