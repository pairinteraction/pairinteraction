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

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

#include "EigenCompat.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

template <bool is_loading>
struct serializer;

template <>
struct serializer<false> {
    void serialize(cereal::BinaryOutputArchive &ar, const void *data, size_t size,
                   const char * /* name */) {
        ar.saveBinary(data, size);
    }
    void serialize(cereal::JSONOutputArchive &ar, const void *data, size_t size, const char *name) {
        ar.saveBinaryValue(data, size, name);
    }
};

template <>
struct serializer<true> {
    void serialize(cereal::BinaryInputArchive &ar, void *data, size_t size,
                   const char * /* name */) {
        ar.loadBinary(data, size);
    }
    void serialize(cereal::JSONInputArchive &ar, void *data, size_t size, const char *name) {
        ar.loadBinaryValue(data, size, name);
    }
};

namespace cereal {

template <class Archive, typename _Scalar, int _Options, typename _Index>
void serialize(Archive &ar, Eigen::SparseMatrix<_Scalar, _Options, _Index> &m,
               unsigned int /* version */) {
    _Index innerSize;
    _Index outerSize;
    _Index valuesSize;

    if (Archive::is_saving::value) {
        innerSize = m.innerSize();
        outerSize = m.outerSize();
        valuesSize = m.nonZeros();

        m.makeCompressed();
    }

    ar &CEREAL_NVP(innerSize);
    ar &CEREAL_NVP(outerSize);
    ar &CEREAL_NVP(valuesSize);

    if (Archive::is_loading::value) {
        _Index rows = (m.IsRowMajor) ? outerSize : innerSize;
        _Index cols = (m.IsRowMajor) ? innerSize : outerSize;
        m.resize(rows, cols);
        m.resizeNonZeros(valuesSize);
    }

    serializer<Archive::is_loading::value> s;
    s.serialize(ar, m.innerIndexPtr(), valuesSize * sizeof(_Index), "innerIndexPtr");
    s.serialize(ar, m.outerIndexPtr(), (outerSize + 1) * sizeof(_Index), "outerIndexPtr");
    s.serialize(ar, m.valuePtr(), valuesSize * sizeof(_Scalar), "valuePtr");

    if (Archive::is_loading::value) {
        m.finalize();
    }
}

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows,
          int _MaxCols>
void serialize(Archive &ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> &m,
               unsigned int /* version */) {
    int rows;
    int cols;

    if (Archive::is_saving::value) {
        rows = m.rows();
        cols = m.cols();
    }

    ar &CEREAL_NVP(rows);
    ar &CEREAL_NVP(cols);

    if (Archive::is_loading::value) {
        m.resize(rows, cols);
    }

    serializer<Archive::is_loading::value> s;
    s.serialize(ar, m.data(), rows * cols * sizeof(_Scalar), "data");
}

} // namespace cereal

#endif // SERIALIZATION_EIGEN_H
