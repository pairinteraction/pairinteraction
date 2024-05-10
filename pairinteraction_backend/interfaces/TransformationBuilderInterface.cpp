#include "interfaces/TransformationBuilderInterface.hpp"
#include "enums/TransformationType.hpp"
#include "utils/euler.hpp"

template <typename Scalar>
Transformation<Scalar>
TransformationBuilderInterface<Scalar>::get_rotator(std::array<real_t, 3> to_z_axis,
                                                    std::array<real_t, 3> to_y_axis) const {
    auto euler_zyz_angles = euler::get_euler_angles(to_z_axis, to_y_axis);
    return this->get_rotator(euler_zyz_angles[0], euler_zyz_angles[1], euler_zyz_angles[2]);
}

// Explicit instantiations
template class TransformationBuilderInterface<float>;
template class TransformationBuilderInterface<double>;
template class TransformationBuilderInterface<std::complex<float>>;
template class TransformationBuilderInterface<std::complex<double>>;
