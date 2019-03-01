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

double GreenTensor::KDelta(int i, int j) {
    if (i == j) {
        return 1.;
    } else {
        return 0.;
    }
}

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
    
    Eigen::Matrix<double, 3, 3> plate_tensor;
    Eigen::Matrix<double, 3, 3> plate_tensor2;
    plate_tensor1 <<  -2.*x*x+zp*zp, 0, -3.*x*zp,
                    0, x*x+zp*zp, 0,
                    3.*x*zp, 0, -x*x+2*zp*zp;    
    dd_tensor += plate_tensor*std::pow(rp,-5.);
}

// def freiraummatrixdq(rho):
//     return np.array((-3./rho**4)*np.tensordot(np.eye(3),irv,axes=0) + (9./rho**4)*np.tensordot(np.tensordot(irv,irv,axes=0),irv,axes=0) - (3./rho**3)*np.tensordot(irv,Amatrix,axes=0) - (3./rho**3)*matrixAikrhoj,dtype=np.complex64)

void GreenTensor::vacuumDipoleQuadrupole(double x, double y, double z) {//TODO check sign according to sign of dipole-dipole tensors.
    Eigen::Matrix<double, 3, 3> Eye = Eigen::Matrix<double, 3, 3>::Identity(3, 3);
    Eigen::Matrix<double, 3, 1> distance;
    distance << x, y, z;
    Eigen::Matrix<double, 3, 3> Amatrix;
    Amatrix << y* y + z * z, -x*y, -x * z, -x*y, x*x+z*z ,-y*z, -x * z, -y*z, y*y + z * z;
    Amatrix = Amatrix/std::pow(distance.norm(),3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                qd_tensor(i, j, k) = 3. / std::pow(distance.norm(), 4.) *
                        (distance(i) *Eye(j,k)  - 3. * distance(i) * distance(j) * distance(k)) +
                    3. / std::pow(distance.norm(), 3.) *
                        (Amatrix(i, j) * distance(k) + Amatrix(i, k) * distance(j));

                dq_tensor(i, j, k) = -3. / std::pow(distance, 4.) *
                        (distance(i) * Eye(j, k) + 3. * distance(i) * distance(j) * distance(k)) +
                    3. / std::pow(distance.norm(), 3.) *
                        (Amatrix(i, j) * distance(k) + Amatrix(i, k) * distance(j));
            }
        }
    }
}

void GreenTensor::plateDipoleQuadrupole(double x, double zA, double zB) {//TODO check sign according to sign of dipole-dipole tensors.
    double zp = zA + zB;
    double rp = std::sqrt(x * x + zp * zp);
    // add something to Eigen/Tensor via << operator?
    // TODO Check for additional 4*pi -> necessary or not? TODO Check sign!
    qd_tensor(0, 0, 0) -= (-6. * x * x * x + 9. * x * zp * zp) / std::pow(rp, 7.);
    qd_tensor(0, 0, 2) -= (12. * x * x * zp - 3. * zp * zp * zp) / std::pow(rp, 7.);
    qd_tensor(0, 1, 1) -= (3. * x * x * x + 3. * x * zp * zp) / std::pow(rp, 7.);
    qd_tensor(0, 2, 0) -= (-12. * x * x * zp + 3. * zp * zp * zp) / std::pow(rp, 7.);
    qd_tensor(0, 2, 2) -= (-3. * x * x * x + 12. * x * zp * zp) / std::pow(rp, 7.);
    qd_tensor(2, 0, 0) -= (-12. * x * x * zp + 3. * zp * zp * zp) / std::pow(rp, 7.);
    qd_tensor(2, 0, 2) -= (-3. * x * x * x + 12. * x * zp * zp) / std::pow(rp, 7.);
    qd_tensor(2, 1, 1) -= (3. * zp * rp * rp) / std::pow(rp, 7.);
    qd_tensor(2, 2, 0) -= (3. * x * x * x - 12. * x * zp * zp) / std::pow(rp, 7.);
    qd_tensor(2, 2, 2) -= (-9. * x * x * zp + zp * zp * zp) / std::pow(rp, 7.);

    dq_tensor(0, 0, 0) -= (6. * x * x * x - 9. * x * zp * zp) / std::pow(rp, 7.);
    dq_tensor(0, 0, 2) -= (-12. * x * x * zp + 3. * zp * zp * zp) / std::pow(rp, 7.);
    dq_tensor(0, 2, 0) -= (12. * x * x * zp - 3. * zp * zp * zp) / std::pow(rp, 7.);
    dq_tensor(0, 2, 2) -= (-3. * x * x * x + 12. * x * zp * zp) / std::pow(rp, 7.);
    dq_tensor(1, 1, 0) -= (-3. * x * x * x + -3. * x * zp * zp) / std::pow(rp, 7.);
    dq_tensor(1, 1, 2) -= (3. * zp * rp * rp) / std::pow(rp, 7.);
    dq_tensor(2, 0, 0) -= (12. * x * x * zp - 3. * zp * zp * zp) / std::pow(rp, 7.);
    dq_tensor(2, 0, 2) -= (3. * x * x * x - 12. * x * zp * zp) /
        std::pow(rp,7.); // TODO check this element, probably wrong in old code (3*x*x - 12 *x*zp*zp)
    dq_tensor(2, 2, 0) -= (-3. * x * x * x + 12 * x * zp * zp) / std::pow(rp, 7.);
    dq_tensor(2, 2, 2) -= (-9. * x * zp * zp + zp * zp * zp) / std::pow(rp, 7.);
}



const Eigen::Matrix<double, 3, 3> &GreenTensor::getDDTensor() {
    if (!dd_tensor_calculated) {

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dd_tensor(i, j) = 0.; // TODO use function from Eigen
            }
        }
        if (surface_distance != std::numeric_limits<double>::max()) {
            x = sqrt(x * x + y * y);
            y = 0.;
            vacuum(x,y,z);
            angle = std::atan(x / z);
            zA = surface_distance - z * std::sin(angle) / 2.;
            zB = surface_distance + z * std::sin(angle) / 2.;
            plate(x, zA, zB);
        }
        else{
            vacuum(x,y,z);
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
