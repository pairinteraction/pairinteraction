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
#include <cereal/cereal.hpp>

#include "EigenCompat.hpp"
#include <Eigen/Core>
#include <Eigen/SparseCore>

template <typename T, typename S>
struct ArrayData {
    T *data;
    S size;

    template <typename Archive>
    void CEREAL_SERIALIZE_FUNCTION_NAME(Archive &ar) {
        using namespace cereal;
        ar(make_size_tag(static_cast<size_type>(size)));
        if constexpr (cereal::traits::is_text_archive<Archive>::value) { // NOLINT
            for (S i = 0; i < size; ++i) {
                ar(data[i]);
            }
        } else { // NOLINT
            ar(binary_data(data, size * sizeof(*data)));
        }
    }
};

template <typename T, typename S>
ArrayData(T, S) -> ArrayData<std::remove_pointer_t<T>, std::decay_t<S>>;

namespace cereal {

template <typename Archive, typename Scalar, int Options, typename StorageIndex>
void CEREAL_SERIALIZE_FUNCTION_NAME(Archive &ar,
                                    Eigen::SparseMatrix<Scalar, Options, StorageIndex> &m) {
    StorageIndex innerSize;
    StorageIndex outerSize;
    StorageIndex valuesSize;

    if constexpr (Archive::is_saving::value) { // NOLINT
        m.makeCompressed();
        innerSize = m.innerSize();
        outerSize = m.outerSize();
        valuesSize = m.nonZeros();
    }

    ar(CEREAL_NVP(innerSize));
    ar(CEREAL_NVP(outerSize));
    ar(CEREAL_NVP(valuesSize));

    if constexpr (Archive::is_loading::value) { // NOLINT
        StorageIndex const rows = (m.IsRowMajor) ? outerSize : innerSize;
        StorageIndex const cols = (m.IsRowMajor) ? innerSize : outerSize;
        m.resize(rows, cols);
        m.resizeNonZeros(valuesSize);
    }

    ar(make_nvp("innerIndexPtr", ArrayData{m.innerIndexPtr(), valuesSize}));
    ar(make_nvp("outerIndexPtr", ArrayData{m.outerIndexPtr(), outerSize + 1}));
    ar(make_nvp("valuePtr", ArrayData{m.valuePtr(), valuesSize}));

    if constexpr (Archive::is_loading::value) { // NOLINT
        m.finalize();
    }
}

template <class Archive, typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void CEREAL_SERIALIZE_FUNCTION_NAME(
    Archive &ar, Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> &m) {
    Eigen::Index rows;
    Eigen::Index cols;

    if constexpr (Archive::is_saving::value) { // NOLINT
        rows = m.rows();
        cols = m.cols();
    }

    ar(CEREAL_NVP(rows));
    ar(CEREAL_NVP(cols));

    if constexpr (Archive::is_loading::value) { // NOLINT
        m.resize(rows, cols);
    }

    ar(make_nvp("data", ArrayData{m.data(), rows * cols}));
}

} // namespace cereal

#endif // SERIALIZATION_EIGEN_H
