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

#ifndef GREENTENSOR_H
#define GREENTENSOR_H

#include "dtypes.hpp"

#include <unsupported/Eigen/CXX11/Tensor>

#include <cmath>
#include <complex>

class GreenTensor {

public:
    using TensorType = Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>>;
    GreenTensor(double x, double y, double z);
    void addSurface(double d);
    const Eigen::Matrix3<double> &getDDTensor();
    const TensorType &getDQTensor();
    const TensorType &getQDTensor();

private:
    Eigen::Matrix3<double> getDDTensorVacuum(double x, double y, double z) const;
    TensorType getDQTensorVacuum(double x, double y, double z) const;
    TensorType getQDTensorVacuum(double x, double y, double z) const;

    Eigen::Matrix3<double> getDDTensorPlate(double x, double zA, double zB) const;
    TensorType getDQTensorPlate(double x, double zA, double zB) const;
    TensorType getQDTensorPlate(double x, double zA, double zB) const;

    Eigen::Matrix3<double> dd_tensor;
    TensorType qd_tensor;
    TensorType dq_tensor;

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
