// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "pairinteraction/utils/eigen_assertion.hpp"

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
template <typename Type, int Size>
using Vector = Matrix<Type, Size, 1>;
} // namespace Eigen
#endif
