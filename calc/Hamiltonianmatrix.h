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

#ifndef HAMILTONIANMATRIX_H
#define HAMILTONIANMATRIX_H

#include "dtypes.h"
#include "Serializable.h"
#include "Basisnames.h"
#include "State.h"
#include "utils.h"

#include <cmath>
#include <cstdint>
#include <stdio.h>
#include <memory>
#include <exception>
#include <stdexcept>
#include <vector>

const uint8_t csr_not_csc = 0x01; // xxx0: csc, xxx1: csr
const uint8_t complex_not_real = 0x02; // xx0x: real, xx1x: complex

class Hamiltonianmatrix : public Serializable {
public:
    Hamiltonianmatrix();
    Hamiltonianmatrix(eigen_sparse_t entries, eigen_sparse_t basis);
    Hamiltonianmatrix(size_t szBasis, size_t szEntries);
    eigen_sparse_t& entries();
    const eigen_sparse_t& entries() const;
    eigen_sparse_t& basis();
    const eigen_sparse_t& basis() const;
    size_t num_basisvectors() const;
    size_t num_coordinates() const;
    void addBasis(idx_t row, idx_t col, scalar_t val);
    void addEntries(idx_t row, idx_t col, scalar_t val);
    void compress(size_t nBasis, size_t nCoordinates);
    std::vector<Hamiltonianmatrix> findSubs() const;
    Hamiltonianmatrix abs() const;
    Hamiltonianmatrix changeBasis(eigen_sparse_t basis) const;
    void applyCutoff(double cutoff);
    void findUnnecessaryStates(std::vector<bool> &isNecessaryCoordinate) const;
    void removeUnnecessaryBasisvectors(const std::vector<bool> &isNecessaryCoordinate);
    void removeUnnecessaryBasisvectors();
    void removeUnnecessaryStates(const std::vector<bool> &isNecessaryCoordinate);
    Hamiltonianmatrix getBlock(const std::vector<ptrdiff_t> &indices);
    void diagonalize();
    friend Hamiltonianmatrix operator+(Hamiltonianmatrix lhs, const Hamiltonianmatrix& rhs);
    friend Hamiltonianmatrix operator-(Hamiltonianmatrix lhs, const Hamiltonianmatrix& rhs);
    friend Hamiltonianmatrix operator*(const scalar_t& lhs,  Hamiltonianmatrix rhs);
    friend Hamiltonianmatrix operator*(Hamiltonianmatrix lhs,  const scalar_t& rhs);
    Hamiltonianmatrix& operator+=(const Hamiltonianmatrix& rhs);
    Hamiltonianmatrix& operator-=(const Hamiltonianmatrix& rhs);
    bytes_t& serialize();
    void doSerialization();
    void deserialize(bytes_t &bytesin);
    void doDeserialization();
    uint64_t hashEntries();
    uint64_t hashBasis();
    void save(std::string fname);
    bool load(std::string fname);
    friend Hamiltonianmatrix combine(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const double &deltaE, std::shared_ptr<BasisnamesTwo> basis_two, const Symmetry &sym);
    friend void energycutoff(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const double &deltaE, std::vector<bool> &necessary);

    template<typename T, typename std::enable_if<utils::is_complex<T>::value>::type* = nullptr>
    void mergeComplex(std::vector<storage_double>& real, std::vector<storage_double>& imag, std::vector<T>& complex) {
        std::vector<storage_double>::iterator real_it, imag_it;
        complex.reserve(real.size());
        for (real_it = real.begin(), imag_it = imag.begin(); real_it != real.end(); ++real_it, ++imag_it) {
            complex.push_back(T(*real_it,*imag_it));
        }
    }

    template<typename T, typename std::enable_if<!utils::is_complex<T>::value>::type* = nullptr>
    void mergeComplex(std::vector<storage_double>& real, std::vector<storage_double>& imag, std::vector<T>& complex) {
        (void) imag;
        complex = real;
    }

    template<typename T, typename std::enable_if<utils::is_complex<T>::value>::type* = nullptr>
    void splitComplex(std::vector<storage_double>& real, std::vector<storage_double>& imag, std::vector<T>& complex) {
        real.reserve(complex.size());
        imag.reserve(imag.size());
        for (auto complex_it = complex.begin(); complex_it != complex.end(); ++complex_it) {
            real.push_back(complex_it->real());
            imag.push_back(complex_it->imag());
        }
    }

    template<typename T, typename std::enable_if<!utils::is_complex<T>::value>::type* = nullptr>
    void splitComplex(std::vector<storage_double>& real, std::vector<storage_double>& imag, std::vector<T>& complex) {
        imag = std::vector<storage_double>();
        real = complex; //std::vector<storage_double>(complex.begin(),complex.end());
    }

protected:
    eigen_sparse_t entries_;
    eigen_sparse_t basis_;
    bytes_t bytes;
    std::vector<eigen_triplet_t> triplets_basis;
    std::vector<eigen_triplet_t> triplets_entries;
};

#endif // HAMILTONIANMATRIX_H
