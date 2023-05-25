/*
 * Copyright (c) 2023 Sebastian Weber, Henri Menke. All rights reserved.
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

#pragma once

#include <Eigen/Core>

// Backport template aliases to older Eigen versions
// https://eigen.tuxfamily.org/dox/group__matrixtypedefs.html
#if !EIGEN_VERSION_AT_LEAST(3, 4, 0)
namespace Eigen {
template <typename Type>
using MatrixX = Matrix<Type, Dynamic, Dynamic>;
template <typename Type>
using VectorX = Matrix<Type, Dynamic, 1>;
template <typename Type>
using Matrix3 = Matrix<Type, 3, 3>;
template <typename Type>
using Vector3 = Matrix<Type, 3, 1>;
} // namespace Eigen
#endif
