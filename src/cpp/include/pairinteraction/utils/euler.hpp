// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"
#include "pairinteraction/utils/eigen_compat.hpp"

#include <Eigen/Dense>
#include <array>
#include <type_traits>

namespace pairinteraction::euler {
/**
 * @function get_rotation_matrix
 *
 * @brief Build a matrix that rotates the coordinate system to the new z-axis and y-axis.
 *
 * @tparam Real The type of a real number.
 *
 * @param to_z_axis The new z-axis.
 *
 * @param to_y_axis The new y-axis.
 *
 * @return The passive rotation matrix.
 */

template <typename Real>
inline Eigen::Matrix3<Real> get_rotation_matrix(std::array<Real, 3> to_z_axis,
                                                std::array<Real, 3> to_y_axis) {
    static_assert(std::is_floating_point_v<Real>);

    auto to_z_axis_mapped = Eigen::Map<Eigen::Vector3<Real>>(to_z_axis.data()).normalized();
    auto to_y_axis_mapped = Eigen::Map<Eigen::Vector3<Real>>(to_y_axis.data()).normalized();

    auto scale = to_z_axis_mapped.norm() * to_y_axis_mapped.norm();
    auto numerical_precision = 100 * scale * std::numeric_limits<Real>::epsilon();

    if (std::abs(to_z_axis_mapped.dot(to_y_axis_mapped)) > numerical_precision) {
        throw std::runtime_error("The z-axis and the y-axis are not orhogonal.");
    }

    Eigen::Matrix3<Real> rotator;
    rotator << to_y_axis_mapped.cross(to_z_axis_mapped), to_y_axis_mapped, to_z_axis_mapped;

    return rotator;
}

/**
 * @function get_euler_angles
 *
 * @brief Extract the Euler angles alpha, beta, gamma
 *
 * @tparam Real The type of a real number.
 *
 * @param to_z_axis The new z-axis.
 *
 * @param to_y_axis The new y-axis.
 *
 * @return The Euler angles.
 */

template <typename Real>
inline std::array<Real, 3> get_euler_angles(std::array<Real, 3> to_z_axis,
                                            std::array<Real, 3> to_y_axis) {
    static_assert(std::is_floating_point_v<Real>);

    auto rotator = get_rotation_matrix(to_z_axis, to_y_axis);
    std::array<Real, 3> euler_zyz{};
    Eigen::Map<Eigen::Vector3<Real>>(euler_zyz.data()) = rotator.eulerAngles(2, 1, 2);
    return euler_zyz;
}

} // namespace pairinteraction::euler
