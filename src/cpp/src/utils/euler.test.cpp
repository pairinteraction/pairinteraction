// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "pairinteraction/utils/euler.hpp"

#include <array>
#include <doctest/doctest.h>

namespace pairinteraction {
DOCTEST_TEST_CASE("construction of rotation matrixes") {
    auto rotator = euler::get_rotation_matrix<double>({0, 0, 1}, {0, 1, 0});
    auto rotator_reference = Eigen::Matrix3<double>::Identity();
    DOCTEST_CHECK((rotator - rotator_reference).norm() == 0);

    rotator = euler::get_rotation_matrix<double>({0, 0, 1}, {1, 1, 0});
    auto y_axis = Eigen::Vector3<double>{0, 1, 0};
    auto rotated_y_axis = rotator * y_axis;
    auto rotated_y_axis_reference = Eigen::Vector3<double>{1, 1, 0}.normalized();
    DOCTEST_CHECK((rotated_y_axis - rotated_y_axis_reference).norm() == 0);

    rotator = euler::get_rotation_matrix<double>({1, 0, 0}, {0, 1, 0});
    auto z_axis = Eigen::Vector3<double>{0, 0, 1};
    auto rotated_z_axis = rotator * z_axis;
    auto rotated_z_axis_reference = Eigen::Vector3<double>{1, 0, 0};
    DOCTEST_CHECK((rotated_z_axis - rotated_z_axis_reference).norm() == 0);

    std::string error_msg = "The z-axis and the y-axis are not orhogonal.";
    DOCTEST_CHECK_THROWS_WITH_AS(euler::get_rotation_matrix<double>({0, 0, 1}, {0, 1, 1});
                                 , error_msg.c_str(), std::runtime_error);
}

DOCTEST_TEST_CASE("construction of zyz euler angles") {
    constexpr double PI = 3.141592653589793238462643383279502884;

    auto euler_angles = euler::get_euler_angles<double>({1, 1, 0}, {0, 0, 1});
    std::array<double, 3> euler_angles_reference{0.25 * PI, 0.5 * PI, 0.5 * PI};
    DOCTEST_CHECK(std::abs(euler_angles[0] - euler_angles_reference[0]) +
                      std::abs(euler_angles[1] - euler_angles_reference[1]) +
                      std::abs(euler_angles[2] - euler_angles_reference[2]) <=
                  1e-6);
}
} // namespace pairinteraction
