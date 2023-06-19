/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef SERIALIZATION_EIGEN_H
#define SERIALIZATION_EIGEN_H

#include "EigenCompat.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace boost {
namespace serialization {

template <class Archive, typename _Scalar, int _Options, typename _Index>
void serialize(Archive &ar, Eigen::SparseMatrix<_Scalar, _Options, _Index> &m,
               const unsigned int /*version*/) {

    _Index innerSize;
    _Index outerSize;
    _Index valuesSize;

    if (Archive::is_saving::value) {
        innerSize = m.innerSize();
        outerSize = m.outerSize();
        valuesSize = m.nonZeros();

        m.makeCompressed();
    }

    ar &innerSize;
    ar &outerSize;
    ar &valuesSize;

    if (Archive::is_loading::value) {
        _Index rows = (m.IsRowMajor) ? outerSize : innerSize;
        _Index cols = (m.IsRowMajor) ? innerSize : outerSize;
        m.resize(rows, cols);
        m.resizeNonZeros(valuesSize);
    }

    ar &make_array(m.innerIndexPtr(), valuesSize);
    ar &make_array(m.outerIndexPtr(), outerSize + 1);
    ar &make_array(m.valuePtr(), valuesSize);

    if (Archive::is_loading::value) {
        m.finalize();
    }
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows,
          int _MaxCols>
void serialize(Archive &ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> &m,
               const unsigned int /*version*/) {

    int rows;
    int cols;

    if (Archive::is_saving::value) {
        rows = m.rows();
        cols = m.cols();
    }

    ar &rows;
    ar &cols;

    if (Archive::is_loading::value) {
        m.resize(rows, cols);
    }

    ar &make_array(m.data(), rows * cols);
}

} // namespace serialization
} // namespace boost

#endif // SERIALIZATION_EIGEN_H
