/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke, Johannes Block. All rights reserved.
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
