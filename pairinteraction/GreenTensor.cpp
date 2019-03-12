/*
 * Copyright (c) 2019 Sebastian Weber, Henri Menke, Johannes Block. All rights reserved.
 *
 * This file is part of the pairinteraction library.
 *
 * The pairinteraction library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The pairinteraction library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with the pairinteraction library. If not, see <http://www.gnu.org/licenses/>.
 */

#include "GreenTensor.h"
#include <limits>

GreenTensor::GreenTensor(double x, double y, double z)
    : x(x), y(y), z(z), zA(std::numeric_limits<double>::max()),
      zB(std::numeric_limits<double>::max()), dd_tensor_calculated(false),
      dq_tensor_calculated(false), qd_tensor_calculated(false) {}

void GreenTensor::addSurface(double d) {
    if (y != 0) {
        throw std::runtime_error("The atoms must be in the xz-plane if a surface is present");
    }

    double angle = std::atan(x / z);

    zA = d - z * std::sin(angle) / 2.;
    zB = d + z * std::sin(angle) / 2.;

    if (zA < 0 || zB < 0) {
        throw std::runtime_error(
            "zA or zB < 0. One of the atoms is inside the plate. Plate is half-space z < 0.");
    }

    dd_tensor_calculated = false;
    dq_tensor_calculated = false;
    qd_tensor_calculated = false;
}

const eigen_matrix33 &GreenTensor::getDDTensor() {
    if (!dd_tensor_calculated) {
        dd_tensor = getDDTensorVacuum(x, y, z);
        if (zA != std::numeric_limits<double>::max()) {
            dd_tensor += getDDTensorPlate(x, zA, zB);
        }
        dd_tensor_calculated = true;
    }

    return dd_tensor;
}

const eigen_tensor333 &GreenTensor::getDQTensor() {
    if (!dq_tensor_calculated) {
        dq_tensor = getDQTensorVacuum(x, y, z);
        if (zA != std::numeric_limits<double>::max()) {
            dq_tensor += getDQTensorPlate(x, zA, zB);
        }
        dq_tensor_calculated = true;
    }
    return dq_tensor;
}

const eigen_tensor333 &GreenTensor::getQDTensor() {
    if (!qd_tensor_calculated) {
        qd_tensor = getQDTensorVacuum(x, y, z);
        if (zA != std::numeric_limits<double>::max()) {
            qd_tensor += getQDTensorPlate(x, zA, zB);
        }
        qd_tensor_calculated = true;
    }
    return qd_tensor;
}

eigen_matrix33 GreenTensor::getDDTensorVacuum(double x, double y, double z) const {
    // Build distance vector
    Eigen::Matrix<double, 3, 1> distance;
    distance << x, y, z;

    // Construct Green tensor
    eigen_matrix33 vacuum_tensor =
        -Eigen::MatrixXd::Identity(3, 3) / std::pow(distance.norm(), 3.) +
        3. * distance * distance.transpose() / std::pow(distance.norm(), 5.);

    return vacuum_tensor;
}

eigen_matrix33 GreenTensor::getDDTensorPlate(double x, double zA, double zB) const {
    // Calculate distances to mirror dipole
    double zp = zA + zB;
    double rp = std::sqrt(x * x + zp * zp);

    // Construct Green tensor
    eigen_matrix33 plate_tensor_second_matrix;
    plate_tensor_second_matrix << x * x, 0., -x * zp, 0., 0., 0., x * zp, 0., x * x;

    eigen_matrix33 plate_tensor =
        Eigen::Vector3d({1., 1., 2.}).asDiagonal().toDenseMatrix() / std::pow(rp, 3.) -
        3. * plate_tensor_second_matrix / std::pow(rp, 5);

    return plate_tensor;
}

eigen_tensor333 GreenTensor::getDQTensorVacuum(double x, double y, double z) const {
    eigen_tensor333 vacuum_tensor;
    vacuum_tensor.setZero();

    eigen_matrix33 Eye = eigen_matrix33::Identity(3, 3);
    Eigen::Matrix<double, 3, 1> distance;
    distance << x, y, z;
    double dist = sqrt(x * x + y * y + z * z);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                vacuum_tensor(i, j, k) += -3 / std::pow(dist, 5.) *
                        (Eye(i, j) * distance(k) + distance(i) * Eye(j, k) +
                         distance(j) * Eye(i, k)) +
                    15. / std::pow(dist, 7.) * distance(i) * distance(j) * distance(k);
            }
        }
    }

    return vacuum_tensor;
}

eigen_tensor333 GreenTensor::getQDTensorVacuum(double x, double y, double z) const {
    eigen_tensor333 vacuum_tensor;
    vacuum_tensor.setZero();

    eigen_matrix33 Eye = eigen_matrix33::Identity(3, 3);
    Eigen::Matrix<double, 3, 1> distance;
    distance << x, y, z;
    double dist = sqrt(x * x + y * y + z * z);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {

                vacuum_tensor(i, j, k) += 3 / std::pow(dist, 5.) *
                        (distance(i) * Eye(j, k) + Eye(i, j) * distance(k) +
                         Eye(i, k) * distance(j)) -
                    15. / std::pow(dist, 7.) * distance(i) * distance(j) * distance(k);
            }
        }
    }

    return vacuum_tensor;
}

eigen_tensor333 GreenTensor::getDQTensorPlate(double x, double zA, double zB) const {
    eigen_tensor333 plate_tensor;
    plate_tensor.setZero();

    double zp = zA + zB;
    Eigen::Matrix<double, 3, 1> distanceplus;
    distanceplus << -x, 0., zp;
    double rp = distanceplus.norm();

    eigen_matrix33 plate_tensor_first_matrix =
        Eigen::Vector3d({1., 1., 2.}).asDiagonal().toDenseMatrix();
    eigen_matrix33 plate_tensor_second_matrix;
    plate_tensor_second_matrix << x * x, 0., -x * zp, 0., 0., 0., x * zp, 0., x * x;

    eigen_tensor333 plate_tensor_second_gradient;
    plate_tensor_second_gradient.setZero();
    plate_tensor_second_gradient(0, 0, 0) = -2. * x;
    plate_tensor_second_gradient(2, 2, 0) = -2. * x;
    plate_tensor_second_gradient(0, 2, 0) = zp;
    plate_tensor_second_gradient(2, 0, 0) = -zp;
    plate_tensor_second_gradient(0, 2, 2) = -x;
    plate_tensor_second_gradient(2, 0, 2) = x;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                plate_tensor(i, j, k) +=
                    -3. * plate_tensor_first_matrix(i, j) * distanceplus(k) / std::pow(rp, 5.);
                plate_tensor(i, j, k) +=
                    15. * plate_tensor_second_matrix(i, j) * distanceplus(k) / std::pow(rp, 7.);
                plate_tensor(i, j, k) +=
                    -3. * plate_tensor_second_gradient(i, j, k) / std::pow(rp, 5.);
            }
        }
    }

    return plate_tensor;
}

eigen_tensor333 GreenTensor::getQDTensorPlate(double x, double zA, double zB) const {
    eigen_tensor333 plate_tensor;
    plate_tensor.setZero();

    double zp = zA + zB;
    Eigen::Matrix<double, 3, 1> distanceplus;
    distanceplus << x, 0., zp;
    double rp = distanceplus.norm();

    eigen_matrix33 plate_tensor_first_matrix =
        Eigen::Vector3d({1., 1., 2.}).asDiagonal().toDenseMatrix();
    eigen_matrix33 plate_tensor_second_matrix;
    plate_tensor_second_matrix << x * x, 0., -x * zp, 0., 0., 0., x * zp, 0., x * x;

    eigen_tensor333 gradient_plate_tensor_second;
    gradient_plate_tensor_second.setZero();
    gradient_plate_tensor_second(0, 0, 0) = 2. * x;
    gradient_plate_tensor_second(0, 2, 2) = 2. * x;
    gradient_plate_tensor_second(0, 0, 2) = -zp;
    gradient_plate_tensor_second(0, 2, 0) = zp;
    gradient_plate_tensor_second(2, 0, 2) = -x;
    gradient_plate_tensor_second(2, 2, 0) = x;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                plate_tensor(i, j, k) +=
                    -3. * distanceplus(i) * plate_tensor_first_matrix(j, k) / std::pow(rp, 5.);
                plate_tensor(i, j, k) +=
                    15. * distanceplus(i) * plate_tensor_second_matrix(j, k) / std::pow(rp, 7.);
                plate_tensor(i, j, k) +=
                    -3. * gradient_plate_tensor_second(i, j, k) / std::pow(rp, 5.);
            }
        }
    }

    return plate_tensor;
}
