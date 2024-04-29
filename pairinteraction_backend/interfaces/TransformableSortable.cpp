#include "interfaces/TransformableSortable.hpp"
#include "utils/euler.hpp"

template <typename Scalar>
Eigen::SparseMatrix<Scalar> TransformableSortable<Scalar>::get_rotator(real_t alpha, real_t beta,
                                                                       real_t gamma) const {
    return this->impl_get_rotator(alpha, beta, gamma);
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar>
TransformableSortable<Scalar>::get_rotator(std::array<real_t, 3> to_z_axis,
                                           std::array<real_t, 3> to_y_axis) const {
    auto euler_zyz_angles = euler::get_euler_angles(to_z_axis, to_y_axis);
    return this->impl_get_rotator(euler_zyz_angles[0], euler_zyz_angles[1], euler_zyz_angles[2]);
}

template <typename Scalar>
std::vector<int> TransformableSortable<Scalar>::get_sorter(SortBy label) const {
    return this->impl_get_sorter(label);
}

template <typename Scalar>
std::vector<int> TransformableSortable<Scalar>::get_blocks(SortBy label) const {
    if (label != sorting) {
        throw std::invalid_argument("The object is not sorted by the requested label.");
    }
    return this->impl_get_blocks(label);
}

template <typename Scalar>
void TransformableSortable<Scalar>::transform(const Eigen::SparseMatrix<Scalar> &transformator) {
    float tolerance = 1e-16;

    Eigen::SparseMatrix<Scalar> identity(transformator.cols(), transformator.cols());
    identity.setIdentity();
    if ((transformator.adjoint() * transformator - identity).norm() > tolerance) {
        throw std::runtime_error("The transformation is not unitary.");
    }

    this->impl_transform(transformator);

    sorting = SortBy::NONE;
    transformation = TransformBy::ARBITRARY;
}

template <typename Scalar>
void TransformableSortable<Scalar>::rotate(real_t alpha, real_t beta, real_t gamma) {
    if (transformation != TransformBy::IDENTITY) {
        throw std::runtime_error("If the object was transformed, it can no longer be rotated.");
    }
    auto rotator = this->impl_get_rotator(alpha, beta, gamma);
    this->impl_transform(rotator);

    sorting = SortBy::NONE;
    transformation = TransformBy::ROTATION;
}

template <typename Scalar>
void TransformableSortable<Scalar>::rotate(std::array<real_t, 3> to_z_axis,
                                           std::array<real_t, 3> to_y_axis) {
    if (transformation != TransformBy::IDENTITY) {
        throw std::runtime_error("If the object was transformed, it can no longer be rotated.");
    }
    auto euler_zyz_angles = euler::get_euler_angles(to_z_axis, to_y_axis);
    auto rotator =
        this->impl_get_rotator(euler_zyz_angles[0], euler_zyz_angles[1], euler_zyz_angles[2]);
    this->impl_transform(rotator);

    sorting = SortBy::NONE;
    transformation = TransformBy::ROTATION;
}

template <typename Scalar>
void TransformableSortable<Scalar>::sort(const std::vector<int> &sorter) {
    this->impl_sort(sorter);

    sorting = SortBy::NONE;
}

template <typename Scalar>
void TransformableSortable<Scalar>::sort(SortBy label) {
    auto sorter = this->impl_get_sorter(label);
    this->impl_sort(sorter);

    sorting = label;
}

// Explicit instantiations
template class TransformableSortable<float>;
template class TransformableSortable<double>;
template class TransformableSortable<std::complex<float>>;
template class TransformableSortable<std::complex<double>>;
