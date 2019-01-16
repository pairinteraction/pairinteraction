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

double GreenTensor::KDelta(int i, int j){
    if(i==j){
        return 1.;
    } else{
        return 0.;
    }
}


void GreenTensor::vacuum(double x, double z) {
    double distance = std::sqrt(x * x + z * z);
    double vecrho[3];
//     z = std::abs(z);
    vecrho[0] = x / (distance);
    vecrho[1] = 0.;
    vecrho[2] = z / (distance);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dd_tensor(i, j) = (KDelta(i,j) - 3. * vecrho[i] * vecrho[j]) / std::pow(distance, 3.);
        }
    }
}

void GreenTensor::vacuumDipoleQuadrupole(double x, double z){
    (void)x;
    (void)z;
}

void GreenTensor::vacuumQuadrupoleDipole(double x, double z){
    (void)x;
    (void)z;
}



// void GreenTensor::vacuumQuadrupoleDipole(double x, double z){
//     double distance = std::sqrt(x * x + z * z);
//     double vecrho[3];
//     vecrho[0] = x / (distance);
//     vecrho[1] = 0.;
//     vecrho[2] = z / (distance);
//     Eigen::Matrix<double, 3, 3> Amatrix;
//     Amatrix << z*z, 0, -x*z,
//                0, distance*distance, 0,
//                -x*z, 0, z*z;
//     Amatrix /= distance*distance;
//     for(int i = 0; i<3; i++){
//         for(int j = 0; j<3; j++){
//             for(int k = 0; k<3; k++){
//                 qd_tensor(i,j,k) = 3./std::pow(distance, 4.) * vecrho[i]*KDelta(j,k) - 9./std::pow(distance,4.)*vecrho[i]*vecrho[j]*vecrho[k] + 3./std::pow(distance,3.)*(Amatrix(i,j)*vecrho[k]  + Amatrix(i,k)*vecrho[j]) ;
//             }
//         }
//     }
// }

// Amatrix = np.array([[ry*ry+rz*rz,-rx*ry,-rx*rz],[-rx*ry,rx*rx+rz*rz,-ry*rz],[-rx*rz,-ry*rz,rx*rx+ry*ry]])/(rx*rx+ry*ry+rz*rz)**(3/2)
// Amatrix = np.array([[rz*rz,0,-rx*rz],[0,rx*rx+rz*rz,0],[-rx*rz,0,rx*rx]])/distance**3

//old Python code:
// irv = identityrvector -> normierter r vector.
// def freiraummatrixqd(rho):
//     return np.array((3./rho**4)*np.tensordot(irv,np.eye(3),axes=0) - (9./rho**4)*np.tensordot(np.tensordot(irv,irv,axes=0),irv,axes=0) + (3./rho**3)*np.tensordot(Amatrix,irv,axes=0) + (3./rho**3)*matrixAikrhoj,dtype=np.complex64)

void GreenTensor::plate(double x, double zA, double zB) {
    if (zA <= 0. || zB <= 0.) {
        throw std::runtime_error(
                "zA or zB <= 0. One of the atoms is inside the plate. Plate is half-space z <= 0.");
    }
    double zp = zA + zB;
    double rp = std::sqrt(x * x + zp * zp);

    dd_tensor(0, 0) -= (1. - (3. * x * x) / (rp * rp)) * std::pow(rp, -3.);
    dd_tensor(0, 2) -= ((3. * x * zp) / (rp * rp)) * std::pow(rp, -3.);
    dd_tensor(1, 1) -= 1. * std::pow(rp, -3.);
    dd_tensor(2, 0) -= (-3. * x * zp / (rp * rp)) * std::pow(rp, -3.);
    dd_tensor(2, 2) -= (2. - 3. * x * x / (rp * rp)) * std::pow(rp, -3.);

}

const Eigen::Matrix<double, 3, 3> &GreenTensor::getDDTensor(){
    if(!dd_tensor_calculated){
    
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dd_tensor(i, j) = 0.; //TODO use function from Eigen
            }
        }
        vacuum(x, z);
        if(surface_distance!=std::numeric_limits<double>::max()){
            zA = surface_distance - z*std::sin(angle)/2.;
            zB = surface_distance + z*std::sin(angle)/2.;
            plate(x,zA,zB);
        }
        
    dd_tensor_calculated = true;
        
    }
    
    return dd_tensor;
}

void GreenTensor::addSurface(double d){
    surface_distance = d;
    dd_tensor_calculated = false;
    dq_tensor_calculated = false;
}
GreenTensor::GreenTensor(double a, double b): x(a), z(b), angle(std::atan(x/z)), zA(std::numeric_limits<double>::max()), zB(std::numeric_limits<double>::max()), dd_tensor_calculated(false), dq_tensor_calculated(false), surface_distance(std::numeric_limits<double>::max()) {    
}


