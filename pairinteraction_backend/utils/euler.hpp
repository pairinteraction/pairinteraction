#pragma once

#include <Eigen/Dense>
#include <array>

namespace euler {

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
Eigen::Matrix<Real, 3, 3> get_rotation_matrix(std::array<Real, 3> to_z_axis,
                                              std::array<Real, 3> to_y_axis) {
    auto to_z_axis_mapped = Eigen::Map<Eigen::Matrix<Real, 3, 1>>(to_z_axis.data()).normalized();
    auto to_y_axis_mapped = Eigen::Map<Eigen::Matrix<Real, 3, 1>>(to_y_axis.data()).normalized();

    Real tolerance = 1e-16;
    if (std::abs(to_z_axis_mapped.dot(to_y_axis_mapped)) > tolerance) {
        throw std::runtime_error("The z-axis and the y-axis are not orhogonal.");
    }

    Eigen::Matrix<Real, 3, 3> rotator;
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
std::array<Real, 3> get_euler_angles(std::array<Real, 3> to_z_axis, std::array<Real, 3> to_y_axis) {
    auto rotator = get_rotation_matrix(to_z_axis, to_y_axis);
    std::array<Real, 3> euler_zyz;
    Eigen::Map<Eigen::Matrix<Real, 3, 1>>(euler_zyz.data()) = rotator.eulerAngles(2, 1, 2);
    return euler_zyz;
}

} // namespace euler

///////////////////////////////////////////////////////////////////////////////////////
// Test cases
///////////////////////////////////////////////////////////////////////////////////////

#include <doctest/doctest.h>
#include <fmt/ostream.h>
#include <spdlog/spdlog.h>

#if FMT_VERSION < 90000
namespace fmt {
template <typename T>
inline auto streamed(T&& v) { return std::forward<T>(v); }
}
#endif

DOCTEST_TEST_CASE("construction of rotation matrixes") {
    auto rotator = euler::get_rotation_matrix<double>({0, 0, 1}, {0, 1, 0});
    auto rotator_reference = Eigen::Matrix<double, 3, 3>::Identity();
    DOCTEST_CHECK((rotator - rotator_reference).norm() == 0);

    rotator = euler::get_rotation_matrix<double>({0, 0, 1}, {1, 1, 0});
    auto y_axis = Eigen::Matrix<double, 3, 1>{0, 1, 0};
    auto rotated_y_axis = rotator * y_axis;
    auto rotated_y_axis_reference = Eigen::Matrix<double, 3, 1>{1, 1, 0}.normalized();
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Rotation matrix:\n{}", fmt::streamed(rotator));
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Rotated y-axis:\n{}", fmt::streamed(rotated_y_axis));
    DOCTEST_CHECK((rotated_y_axis - rotated_y_axis_reference).norm() == 0);

    rotator = euler::get_rotation_matrix<double>({1, 0, 0}, {0, 1, 0});
    auto z_axis = Eigen::Matrix<double, 3, 1>{0, 0, 1};
    auto rotated_z_axis = rotator * z_axis;
    auto rotated_z_axis_reference = Eigen::Matrix<double, 3, 1>{1, 0, 0};
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Rotation matrix:\n{}", fmt::streamed(rotator));
    SPDLOG_LOGGER_INFO(spdlog::get("doctest"), "Rotated z-axis:\n{}", fmt::streamed(rotated_z_axis));
    DOCTEST_CHECK((rotated_z_axis - rotated_z_axis_reference).norm() == 0);

    std::string error_msg = "The z-axis and the y-axis are not orhogonal.";
    DOCTEST_CHECK_THROWS_WITH_AS(euler::get_rotation_matrix<double>({0, 0, 1}, {0, 1, 1});
                                 , error_msg.c_str(), std::runtime_error);
}
