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

#ifndef DTYPES_H
#define DTYPES_H

#include <vector>
#include <array>
#include <Eigen/Sparse>
#include <Eigen/Dense>

typedef std::complex<double> complex_t;
typedef uint32_t idx_t;
typedef double storage_double; // TODO has always to be the same as double
typedef int32_t storage_idx_t;

typedef uint8_t byte_t;
typedef std::vector<byte_t> bytes_t;
typedef std::nullptr_t invalid_t;

#ifdef USE_COMPLEX
typedef complex_t scalar_t;
#else
typedef double scalar_t;
#endif

typedef Eigen::Triplet<scalar_t> eigen_triplet_t;
typedef Eigen::Triplet<double> eigen_triplet_double_t;
typedef Eigen::SparseMatrix<scalar_t> eigen_sparse_t;
typedef Eigen::SparseMatrix<double> eigen_sparse_double_t;
typedef Eigen::SparseMatrix<scalar_t>::InnerIterator eigen_iterator_t;
typedef Eigen::SparseMatrix<double>::InnerIterator eigen_iterator_double_t;
typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,Eigen::Dynamic> eigen_dense_t;
typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,1> eigen_vector_t;
typedef Eigen::Matrix<double,Eigen::Dynamic,1> eigen_vector_double_t;


class Triple {
public:
    Triple() : row(0), col(0), val(0) { }
    Triple(idx_t row, idx_t col, scalar_t val) : row(row), col(col), val(val) { }
    idx_t row;
    idx_t col;
    scalar_t val;
};

enum parity_t {
    NA = INT_MAX,
    EVEN = 1,
    ODD = -1,
};

struct Symmetry {
    parity_t inversion;
    parity_t reflection;
    parity_t permutation;
    int rotation;

    // Comparison operator that is needed if an object of type Symmetry is used as key for std::map
    friend bool operator< (const Symmetry& s1, const Symmetry& s2)
    {
        std::array<int, 5> syms1{ {s1.inversion, s1.reflection, s1.permutation, s1.rotation} };
        std::array<int, 5> syms2{ {s2.inversion, s2.reflection, s2.permutation, s2.rotation} };

        for (size_t i = 0; i < syms1.size(); ++i) {
            if (syms1[i] < syms2[i]) {
                return true;
            } else if (syms1[i] > syms2[i]) {
                return false;
            }
        }
        return false;
    }
};

#endif
