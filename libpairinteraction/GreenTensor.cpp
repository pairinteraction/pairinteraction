/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke, Johannes Block. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "GreenTensor.h"

void GreenTensor::vacuum(double x, double y, double z) {
    // Build distance vector
    Eigen::Matrix<double, 3, 1> distance;
    distance << x, y, z;

    // Construct Green tensor
    dd_tensor = -Eigen::MatrixXd::Identity(3, 3) / std::pow(distance.norm(), 3) +
        3 * distance * distance.transpose() / std::pow(distance.norm(), 5);
}

void GreenTensor::plate(double x, double zA, double zB) {
    if (zA < 0 || zB < 0) {
        throw std::runtime_error(
            "zA or zB < 0. One of the atoms is inside the plate. Plate is half-space z < 0.");
    }

    double zp = zA + zB;
    double rp = std::sqrt(x * x + zp * zp);

    Eigen::Matrix<double, 3, 3> plate_tensor_second_matrix;
    plate_tensor_second_matrix << x * x, 0, -x * zp, 0, 0, 0, x * zp, 0, x * x;

    // Construct Green tensor
    Eigen::Matrix<double, 3, 3> plate_tensor =
        Eigen::Vector3d({1, 1, 2}).asDiagonal().toDenseMatrix() / std::pow(rp, 3) -
        3 * plate_tensor_second_matrix / std::pow(rp, 5);

    // Add Green tensor to the vacuum Green tensor
    dd_tensor += plate_tensor;
}


void GreenTensor::vacuumDipoleQuadrupole(
    double x, double y, double z) {
    Eigen::Matrix<double, 3, 3> Eye = Eigen::Matrix<double, 3, 3>::Identity(3, 3);
    Eigen::Matrix<double, 3, 1> distance;
    distance << x, y, z; //TODO save distance.norm() as r, rho or something. Aim for consistent notation!
    Eigen::Matrix<double, 3, 3> Amatrix;
    Amatrix << y * y + z * z, -x * y, -x * z, -x * y, x * x + z * z, -y * z, -x * z, -y * z,
        y * y + z * z;
    Amatrix = Amatrix / std::pow(distance.norm(), 3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                qd_tensor(i, j, k) = 3. / std::pow(distance.norm(), 4.) *
                        (distance(i)/distance.norm() * Eye(j, k) - 3. * distance(i)/distance.norm()  * distance(j)/distance.norm()  * distance(k)/distance.norm() ) +
                    3. / std::pow(distance.norm(), 3.) *
                        (Amatrix(i, j) * distance(k) + Amatrix(i, k) * distance(j));

                dq_tensor(i, j, k) = 3. / std::pow(distance.norm(), 4.) *
                        ( Eye(i, j)*distance(k) + 3. * distance(i)/distance.norm() * distance(j)/distance.norm() * distance(k)/distance.norm()) -
                    3. / std::pow(distance.norm(), 3.) *
                        (distance(i)/distance.norm() * Amatrix(j,k) + Amatrix(i,k) * distance(j)/distance.norm() );
            }
        }
    }
}

void GreenTensor::plateDipoleQuadrupole(
    double x, double zA, double zB) {
    double zp = zA + zB;
    Eigen::Matrix<double, 3, 1> distanceplus1;
    distanceplus1 << x, 0, zp;
    Eigen::Matrix<double, 3, 1> distanceplus2;
    distanceplus2 << -x, 0, zp;
    double rp = distanceplus1.norm();
    
    Eigen::Matrix<double, 3, 3> plate_tensor_first_matrix = Eigen::Vector3d({1, 1, 2}).asDiagonal().toDenseMatrix();
    Eigen::Matrix<double, 3, 3> plate_tensor_second_matrix;
    plate_tensor_second_matrix << x * x, 0, -x * zp, 0, 0, 0, x * zp, 0, x * x;
    
    Eigen::Tensor<double, 3> gradient_plate_tensor_second;
    gradient_plate_tensor_second(0,0,0) = 2 * x;
    gradient_plate_tensor_second(0,2,2) = 2 * x;
    gradient_plate_tensor_second(0,0,2) = -zp;
    gradient_plate_tensor_second(0,2,0) = zp;
    gradient_plate_tensor_second(2,0,2) = -x;
    gradient_plate_tensor_second(2,2,0) = x;
    
    Eigen::Tensor<double, 3> plate_tensor_second_gradient;
    plate_tensor_second_gradient(0,0,0) = -2 * x;
    plate_tensor_second_gradient(2,2,0) = -2 * x;
    plate_tensor_second_gradient(0,2,0) = zp;
    plate_tensor_second_gradient(2,0,0) = -zp;
    plate_tensor_second_gradient(0,2,2) = -x;
    plate_tensor_second_gradient(2,0,2) = x; 
    
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                qd_tensor(i,j,k) += -3 * distanceplus1(i)*plate_tensor_first_matrix(j,k)/std::pow(rp,5.);
                qd_tensor(i,j,k) += 15 * distanceplus1(i)*plate_tensor_second_matrix(j,k)/std::pow(rp,7.);
                qd_tensor(i,j,k) += -3 * gradient_plate_tensor_second(i,j,k)/std::pow(rp,5.); 
                
                dq_tensor(i,j,k) += -3 * plate_tensor_first_matrix(i,j)*distanceplus2(k)/std::pow(rp,5.);
                dq_tensor(i,j,k) += 15 * plate_tensor_second_matrix(i,j)*distanceplus2(k)/std::pow(rp,7.);
                dq_tensor(i,j,k) += -3 * plate_tensor_second_gradient(i,j,k)/std::pow(rp,5.);
            }
        }
    }
    
}

const Eigen::Matrix<double, 3, 3> &GreenTensor::getDDTensor() {
    if (!dd_tensor_calculated) {
        dd_tensor = Eigen::Matrix<double, 3, 3>::Zero();
        
        if (surface_distance != std::numeric_limits<double>::max()) {
            x = sqrt(x * x + y * y);
            y = 0.;
            vacuum(x, y, z);
            angle = std::atan(x / z);
            zA = surface_distance - z * std::sin(angle) / 2.;
            zB = surface_distance + z * std::sin(angle) / 2.;
            plate(x, zA, zB);
        } else {
            vacuum(x, y, z);
        }

        dd_tensor_calculated = true;
    }

    return dd_tensor;
}

void GreenTensor::addSurface(double d) {
    surface_distance = d;
    dd_tensor_calculated = false;
    dq_tensor_calculated = false;
}
GreenTensor::GreenTensor(double x, double y, double z)
    : x(x), y(y), z(z), angle(std::atan(x / z)), zA(std::numeric_limits<double>::max()),
      zB(std::numeric_limits<double>::max()), dd_tensor_calculated(false),
      dq_tensor_calculated(false), surface_distance(std::numeric_limits<double>::max()) {}
