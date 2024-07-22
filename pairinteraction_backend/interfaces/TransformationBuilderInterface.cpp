#include "interfaces/TransformationBuilderInterface.hpp"

#include "enums/TransformationType.hpp"
#include "utils/euler.hpp"

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

Blocks::Blocks(std::vector<int> start, std::vector<TransformationType> transformation_type)
    : start(std::move(start)), transformation_type(std::move(transformation_type)) {}

Blocks::Blocks(std::vector<int> start)
    : start(std::move(start)), transformation_type({TransformationType::ARBITRARY}) {}

template <typename Scalar>
Transformation<Scalar>
TransformationBuilderInterface<Scalar>::get_rotator(const std::array<real_t, 3> &to_z_axis,
                                                    const std::array<real_t, 3> &to_y_axis) const {
    auto euler_zyz_angles = euler::get_euler_angles(to_z_axis, to_y_axis);
    return this->get_rotator(euler_zyz_angles[0], euler_zyz_angles[1], euler_zyz_angles[2]);
}

// Explicit instantiations
template struct Transformation<float>;
template struct Transformation<double>;
template struct Transformation<std::complex<float>>;
template struct Transformation<std::complex<double>>;

template class TransformationBuilderInterface<float>;
template class TransformationBuilderInterface<double>;
template class TransformationBuilderInterface<std::complex<float>>;
template class TransformationBuilderInterface<std::complex<double>>;
