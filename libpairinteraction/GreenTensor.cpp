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

void GreenTensor::vacuum(double x, double z) {
    double distance = std::sqrt(x * x + z * z);
    double Kdelta;
    double vecrho[3];
    z = std::abs(z);
    vecrho[0] = x / (distance * distance);
    vecrho[1] = 0.;
    vecrho[2] = z / (distance * distance);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Kdelta = 0.;
            if (i == j) {
                Kdelta = 1.;
            }
            tensor(i, j) = (Kdelta - 3. * vecrho[i] * vecrho[j]) / std::pow(distance, 3.);
        }
    }
}

void GreenTensor::plate(double x, double zA, double zB) {
    if (zA < 0. || zB < 0.) {
        std::cout << "error! z<0" << std::endl;
    }
    double zp = zA + zB;
    double rp = std::sqrt(x * x + zp * zp);
    tensor(0, 0) += (1. - (3. * x * x) / (rp * rp)) * std::pow(rp, -3.);
    tensor(0, 2) += ((3. * x * zp) / (rp * rp)) * std::pow(rp, -3.);
    tensor(1, 1) += 1. * std::pow(rp, -3.);
    tensor(2, 0) += (-3. * x * zp / (rp * rp)) * std::pow(rp, -3.);
    tensor(2, 2) += (2. - 3. * x * x / (rp * rp)) * std::pow(rp, -3.);
}
GreenTensor::GreenTensor(double x, double z) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            tensor(i, j) = 0.;
        }
    }
    vacuum(x, z);
}

GreenTensor::GreenTensor(double x, double zA, double zB) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            tensor(i, j) = 0.;
        }
    }
    vacuum(x, zA - zB);
    plate(x, zA, zB);
}
