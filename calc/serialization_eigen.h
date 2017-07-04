/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
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

#ifndef SERIALIZATION_EIGEN_H
#define SERIALIZATION_EIGEN_H

#include <Eigen/Sparse>

namespace boost {
namespace serialization {

template <class Archive, typename _Scalar, int _Options, typename _Index>
void serialize(Archive& ar, Eigen::SparseMatrix<_Scalar,_Options,_Index>& m, const unsigned int version) {
    (void)version;

    _Index innerSize;
    _Index outerSize;
    _Index valuesSize;

    if(Archive::is_saving::value) {
        innerSize = m.innerSize();
        outerSize = m.outerSize();
        valuesSize = m.nonZeros();

        m.makeCompressed();
    }

    ar & innerSize;
    ar & outerSize;
    ar & valuesSize;

    if(Archive::is_loading::value) {
        _Index rows = (m.IsRowMajor)? outerSize : innerSize;
        _Index cols = (m.IsRowMajor)? innerSize : outerSize;
        m.resize(rows, cols);
        m.resizeNonZeros(valuesSize);
    }

    ar & make_array(m.innerIndexPtr(), valuesSize);
    ar & make_array(m.outerIndexPtr(), outerSize);
    ar & make_array(m.valuePtr(), valuesSize);

    if(Archive::is_loading::value) {
        m.finalize();
    }
}

}
}

#endif // SERIALIZATION_EIGEN_H
