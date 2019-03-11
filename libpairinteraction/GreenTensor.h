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

#ifndef GREENTENSOR_H
#define GREENTENSOR_H

#include "dtypes.h"
#include <cmath>
#include <complex>

class GreenTensor {

public:
    GreenTensor(double x, double y, double z);
    void addSurface(double d);
    const eigen_matrix33 &getDDTensor();
    const eigen_tensor333 &getDQTensor();
    const eigen_tensor333 &getQDTensor();

private:
    eigen_matrix33 getDDTensorVacuum(double x, double y, double z) const;
    eigen_tensor333 getDQTensorVacuum(double x, double y, double z) const;
    eigen_tensor333 getQDTensorVacuum(double x, double y, double z) const;

    eigen_matrix33 getDDTensorPlate(double x, double zA, double zB) const;
    eigen_tensor333 getDQTensorPlate(double x, double zA, double zB) const;
    eigen_tensor333 getQDTensorPlate(double x, double zA, double zB) const;

    eigen_matrix33 dd_tensor;
    eigen_tensor333 qd_tensor;
    eigen_tensor333 dq_tensor;

    double x;
    double y;
    double z;
    double zA;
    double zB;

    bool dd_tensor_calculated;
    bool dq_tensor_calculated;
    bool qd_tensor_calculated;
};

#endif
