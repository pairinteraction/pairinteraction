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
#include <iostream>
#include <iomanip>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <mpi.h>

#define REAL_T MPI_FLOAT
#define IDX_T  MPI_UNSIGNED
#define BYTE_T MPI_UNSIGNED_CHAR
constexpr auto ORDER = Eigen::ColMajor;

typedef double real_t; // TODO we should always use double otherwise the diagonalization results are very noisy
typedef std::complex<real_t> complex_t;
typedef uint32_t idx_t;
typedef real_t storage_real_t; // TODO has always to be the same as real_t
typedef int32_t storage_idx_t;
typedef int32_t eigen_idx_t;

typedef uint8_t byte_t;
typedef std::vector<byte_t> bytes_t;
typedef std::nullptr_t invalid_t;

#ifdef USE_COMPLEX
    typedef complex_t scalar_t;
#else
    typedef real_t scalar_t;
#endif

typedef Eigen::Triplet<scalar_t> eigen_triplet_t;
typedef Eigen::Triplet<real_t> eigen_triplet_real_t;
typedef Eigen::SparseMatrix<scalar_t, ORDER, eigen_idx_t> eigen_sparse_t;
typedef Eigen::SparseMatrix<real_t, ORDER, eigen_idx_t> eigen_sparse_real_t;
typedef Eigen::SparseMatrix<scalar_t, ORDER, eigen_idx_t>::InnerIterator eigen_iterator_t;
typedef Eigen::SparseMatrix<real_t, ORDER, eigen_idx_t>::InnerIterator eigen_iterator_real_t;
typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,Eigen::Dynamic,ORDER> eigen_dense_t;
typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,1,ORDER> eigen_vector_t;
typedef Eigen::Matrix<real_t,Eigen::Dynamic,1,ORDER> eigen_vector_real_t;


class Triple {
public:
    Triple() : row(0), col(0), val(0) { }
    Triple(idx_t row, idx_t col, scalar_t val) : row(row), col(col), val(val) { }
    idx_t row;
    idx_t col;
    scalar_t val;
};

#endif
