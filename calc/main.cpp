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

#include "dtypes.h"
#include "Serializable.h"
#include "QuantumDefect.h"
#include "Basisnames.h"
#include "MatrixElements.h"
#include "SQLite.h"
#include "ConfParser.h"
#include "State.h"

#include <memory>
#include <tuple>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <stdio.h>
#include <inttypes.h>
#include <sstream>

#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include <assert.h>

#include <unordered_set>

#include <unordered_map>

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/algorithm/hex.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/program_options.hpp>
#include <boost/functional/hash.hpp>

#include <numeric>

/*
///////////////////// TODOs /////////////////////

* parallelize construction of Hamiltonian
* construct only one half of symmetric matrices
* check why very small Bz fields (e.g. 1e-12) leads to a large basis -> numerical error?

*/

const uint8_t csr_not_csc = 0x01; // xxx0: csc, xxx1: csr
const uint8_t complex_not_real = 0x02; // xx0x: real, xx1x: complex

class Hamiltonianmatrix : public Serializable {
public:
    Hamiltonianmatrix() : Serializable() {}
    //Hamiltonianmatrix(size_t nBasis, size_t nCoordinates) : Serializable(), entries_(nBasis,nBasis), basis_(nCoordinates,nBasis)  {}
    Hamiltonianmatrix(eigen_sparse_t entries, eigen_sparse_t basis) : Serializable(), entries_(entries), basis_(basis) {}

    Hamiltonianmatrix(size_t szBasis, size_t szEntries) : Serializable()  {
        triplets_basis.reserve(szBasis);
        triplets_entries.reserve(szEntries);
    }

    eigen_sparse_t& entries() {
        bytes.clear();
        return entries_;
    }
    const eigen_sparse_t& entries() const {
        return entries_;
    }
    eigen_sparse_t& basis() {
        bytes.clear();
        return basis_;
    }
    const eigen_sparse_t& basis() const {
        return basis_;
    }
    size_t num_basisvectors() const {
        return basis_.cols();
    }
    size_t num_coordinates() const {
        return basis_.rows();
    }

    void addBasis(idx_t row, idx_t col, scalar_t val) {
        triplets_basis.push_back(eigen_triplet_t(row,col,val));
    }
    void addEntries(idx_t row, idx_t col, scalar_t val) {
        triplets_entries.push_back(eigen_triplet_t(row,col,val));
    }
    void compress(size_t nBasis, size_t nCoordinates) {
        basis_.resize(nCoordinates,nBasis);
        entries_.resize(nBasis,nBasis);
        basis_.setFromTriplets(triplets_basis.begin(), triplets_basis.end());
        entries_.setFromTriplets(triplets_entries.begin(), triplets_entries.end());
        triplets_basis.clear();
        triplets_entries.clear();
    }

    std::vector<Hamiltonianmatrix> findSubs() const { // TODO
        std::vector<Hamiltonianmatrix> submatrices;
        submatrices.push_back(*this);
        return submatrices;
    }
    Hamiltonianmatrix abs() const {
        return Hamiltonianmatrix(entries_.cwiseAbs().cast<scalar_t>(), basis_);
    }
    Hamiltonianmatrix changeBasis(eigen_sparse_t basis) const{
        auto transformator = basis_.adjoint()*basis;
        auto entries = transformator.adjoint()*entries_*transformator;
        return Hamiltonianmatrix(entries, basis);
    }

    void applyCutoff(real_t cutoff) {
        bytes.clear();

        // build transformator
        eigen_vector_t diag = entries_.diagonal();

        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(num_basisvectors());

        size_t idxBasis = 0;
        for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
            if (std::abs(diag[idx]) < cutoff) {
                triplets_transformator.push_back(eigen_triplet_t(idx,idxBasis++,1));
            }
        }

        eigen_sparse_t transformator(this->num_basisvectors(),idxBasis);
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // apply transformator
        basis_ = basis_*transformator;
        entries_= transformator.adjoint()*entries_*transformator;
    }

    void findUnnecessaryStates(std::vector<bool> &isNecessaryCoordinate) const {
        std::vector<real_t> isNecessaryCoordinate_real(num_coordinates(),0);
        for (eigen_idx_t k=0; k<basis_.outerSize(); ++k) {
            for (eigen_iterator_t triple(basis_,k); triple; ++triple) {
                isNecessaryCoordinate_real[triple.row()] += std::pow(std::abs(triple.value()),2);
            }
        }

        for (size_t idx = 0; idx < this->num_coordinates(); ++idx) {
            if (isNecessaryCoordinate_real[idx] > 0.05) { // TODO
                isNecessaryCoordinate[idx] = true;
            }
        }
    }

    void removeUnnecessaryBasisvectors(const std::vector<bool> &isNecessaryCoordinate) {
        bytes.clear();

        // build transformator
        std::vector<real_t> isNecessaryBasisvector(num_basisvectors(),0);
        for (eigen_idx_t k_1=0; k_1<basis_.outerSize(); ++k_1) {
            for (eigen_iterator_t triple(basis_,k_1); triple; ++triple) {
                ptrdiff_t col = triple.col(); // basis vector
                ptrdiff_t row = triple.row(); // coordinate
                if (isNecessaryCoordinate[row]) {
                    isNecessaryBasisvector[col] += std::pow(std::abs(triple.value()),2);
                }
            }
        }

        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(num_basisvectors());

        size_t idxBasis = 0;
        for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
            if (isNecessaryBasisvector[idx] > 0.05) { // TODO
                triplets_transformator.push_back(eigen_triplet_t(idx,idxBasis++,1));
            }
        }

        eigen_sparse_t transformator(this->num_basisvectors(),idxBasis);
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // apply transformator
        basis_ = basis_*transformator;
        entries_= transformator.adjoint()*entries_*transformator;
    }

    void removeUnnecessaryBasisvectors() {
        bytes.clear();

        // build transformator
        std::vector<real_t> isNecessaryBasisvector(num_basisvectors(),0);
        for (eigen_idx_t k_1=0; k_1<basis_.outerSize(); ++k_1) {
            for (eigen_iterator_t triple(basis_,k_1); triple; ++triple) {
                ptrdiff_t col = triple.col(); // basis vector
                isNecessaryBasisvector[col] += std::pow(std::abs(triple.value()),2);
            }
        }

        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(num_basisvectors());

        size_t idxBasis = 0;
        for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
            if (isNecessaryBasisvector[idx] > 0.05) {
                triplets_transformator.push_back(eigen_triplet_t(idx,idxBasis++,1));
            }
        }

        eigen_sparse_t transformator(this->num_basisvectors(),idxBasis);
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // apply transformator
        basis_ = basis_*transformator;
        entries_= transformator.adjoint()*entries_*transformator;
    }

    void removeUnnecessaryStates(const std::vector<bool> &isNecessaryCoordinate) {
        bytes.clear();

        // build transformator
        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(num_coordinates());

        size_t idxCoordinate = 0;
        for (size_t idx = 0; idx < this->num_coordinates(); ++idx) {
            if (isNecessaryCoordinate[idx]) {
                triplets_transformator.push_back(eigen_triplet_t(idxCoordinate++,idx,1));
            }
        }

        eigen_sparse_t transformator(idxCoordinate,this->num_coordinates());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // apply transformator
        basis_ = transformator*basis_;
    }

    Hamiltonianmatrix getBlock(const std::vector<ptrdiff_t> &indices) {
        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(indices.size());
        for (size_t idx = 0; idx < indices.size(); ++idx) {
            triplets_transformator.push_back(eigen_triplet_t(indices[idx],idx,1));
        }
        eigen_sparse_t transformator(this->num_basisvectors(),indices.size());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        eigen_sparse_t block_entries = transformator.adjoint()*entries_*transformator;
        eigen_sparse_t block_basis = basis_*transformator;

        return Hamiltonianmatrix(block_entries,  block_basis);
    }

    void diagonalize() {
        if (this->num_basisvectors() > 1) {
            // diagonalization
            Eigen::SelfAdjointEigenSolver<eigen_dense_t> eigensolver(eigen_dense_t(this->entries()));

            // eigenvalues and eigenvectors
            eigen_vector_real_t evals = eigensolver.eigenvalues();
            eigen_sparse_t evecs = eigensolver.eigenvectors().sparseView(1e-4,0.5);

            this->entries().setZero();
            this->entries().reserve(evals.size());
            for (eigen_idx_t idx = 0; idx < evals.size(); ++idx) {
                this->entries().insert(idx, idx) = evals.coeffRef(idx);
            }
            this->entries().makeCompressed();

            this->basis() = (this->basis() * evecs).pruned(1e-4,0.5);
        }
    }

    friend Hamiltonianmatrix operator+(Hamiltonianmatrix lhs, const Hamiltonianmatrix& rhs) {
        lhs.bytes.clear();
        lhs.entries_ += rhs.entries_;
        return lhs;
    }
    friend Hamiltonianmatrix operator-(Hamiltonianmatrix lhs, const Hamiltonianmatrix& rhs) {
        lhs.bytes.clear();
        lhs.entries_ -= rhs.entries_;
        return lhs;
    }
    friend Hamiltonianmatrix operator*(const scalar_t& lhs,  Hamiltonianmatrix rhs) {
        rhs.bytes.clear();
        rhs.entries_ *= lhs;
        return rhs;
    }
    friend Hamiltonianmatrix operator*(Hamiltonianmatrix lhs,  const scalar_t& rhs) {
        lhs.bytes.clear();
        lhs.entries_ *= rhs;
        return lhs;
    }
    Hamiltonianmatrix& operator+=(const Hamiltonianmatrix& rhs) {
        bytes.clear();
        entries_ += rhs.entries_;
        return *this;
    }
    Hamiltonianmatrix& operator-=(const Hamiltonianmatrix& rhs) {
        bytes.clear();
        entries_ -= rhs.entries_;
        return *this;
    }

    bytes_t& serialize() {
        doSerialization();
        return bytes;
    }

    template<typename T, typename std::enable_if<utils::is_complex<T>::value>::type* = nullptr>
    void mergeComplex(std::vector<storage_real_t>& real, std::vector<storage_real_t>& imag, std::vector<T>& complex) {
        std::vector<storage_real_t>::iterator real_it, imag_it;
        complex.reserve(real.size());
        for (real_it = real.begin(), imag_it = imag.begin(); real_it != real.end(); ++real_it, ++imag_it) {
            complex.push_back(T(*real_it,*imag_it));
        }
    }

    template<typename T, typename std::enable_if<!utils::is_complex<T>::value>::type* = nullptr>
    void mergeComplex(std::vector<storage_real_t>& real, std::vector<storage_real_t>& imag, std::vector<T>& complex) {
        (void) imag;
        complex = real;
    }

    template<typename T, typename std::enable_if<utils::is_complex<T>::value>::type* = nullptr>
    void splitComplex(std::vector<storage_real_t>& real, std::vector<storage_real_t>& imag, std::vector<T>& complex) {
        real.reserve(complex.size());
        imag.reserve(imag.size());
        for (auto complex_it = complex.begin(); complex_it != complex.end(); ++complex_it) {
            real.push_back(complex_it->real());
            imag.push_back(complex_it->imag());
        }
    }

    template<typename T, typename std::enable_if<!utils::is_complex<T>::value>::type* = nullptr>
    void splitComplex(std::vector<storage_real_t>& real, std::vector<storage_real_t>& imag, std::vector<T>& complex) {
        imag = std::vector<storage_real_t>();
        real = complex; //std::vector<storage_real_t>(complex.begin(),complex.end());
    }

    void doSerialization() {
        if (bytes.size() == 0) {
            entries_.makeCompressed();
            basis_.makeCompressed();

            // convert matrix "entries" to vectors of primitive data types
            byte_t entries_flags = 0;
            if (entries_.IsRowMajor) {
                entries_flags |= csr_not_csc;
            }
            if (utils::is_complex<scalar_t>::value) {
                entries_flags |= complex_not_real;
            }
            storage_idx_t entries_rows = entries_.rows();
            storage_idx_t entries_cols = entries_.cols();
            std::vector<scalar_t> entries_data(entries_.valuePtr(), entries_.valuePtr()+entries_.nonZeros());
            std::vector<storage_real_t> entries_data_real, entries_data_imag;
            splitComplex(entries_data_real,entries_data_imag,entries_data);
            std::vector<storage_idx_t> entries_indices(entries_.innerIndexPtr(), entries_.innerIndexPtr()+entries_.nonZeros());
            std::vector<storage_idx_t> entries_indptr(entries_.outerIndexPtr(), entries_.outerIndexPtr()+entries_.outerSize());

            // convert matrix "basis" to vectors of primitive data types
            byte_t basis_flags = 0;
            if (basis_.IsRowMajor) {
                basis_flags |= csr_not_csc;
            }
            if (utils::is_complex<scalar_t>::value) {
                basis_flags |= complex_not_real;
            }
            storage_idx_t basis_rows = basis_.rows();
            storage_idx_t basis_cols = basis_.cols();
            std::vector<scalar_t> basis_data(basis_.valuePtr(), basis_.valuePtr()+basis_.nonZeros());
            std::vector<storage_real_t> basis_data_real, basis_data_imag;
            splitComplex(basis_data_real,basis_data_imag,basis_data);
            std::vector<storage_idx_t> basis_indices(basis_.innerIndexPtr(), basis_.innerIndexPtr()+basis_.nonZeros());
            std::vector<storage_idx_t> basis_indptr(basis_.outerIndexPtr(), basis_.outerIndexPtr()+basis_.outerSize());

            // serialize vectors of primitive data types
            Serializer s;
            s << entries_flags;
            s << entries_rows;
            s << entries_cols;
            s << entries_data_real;
            if (entries_flags & complex_not_real) s << entries_data_imag;
            s << entries_indices;
            s << entries_indptr;
            s << basis_flags;
            s << basis_rows;
            s << basis_cols;
            s << basis_data_real;
            if (basis_flags & complex_not_real) s << basis_data_imag;
            s << basis_indices;
            s << basis_indptr;
            s.save(bytes);
        }
    }

    void deserialize(bytes_t &bytesin) {
        bytes = bytesin;
        doDeserialization();
    }

    void doDeserialization() {
        // deserialize vectors of primitive data types
        byte_t entries_flags;
        storage_idx_t entries_rows, entries_cols;
        std::vector<storage_real_t> entries_data_real, entries_data_imag;
        std::vector<idx_t> entries_indices;
        std::vector<idx_t> entries_indptr;
        byte_t basis_flags;
        storage_idx_t basis_rows, basis_cols;
        std::vector<storage_real_t> basis_data_real, basis_data_imag;
        std::vector<idx_t> basis_indices;
        std::vector<idx_t> basis_indptr;

        Serializer s;
        s.load(bytes);
        s >> entries_flags;
        s >> entries_rows;
        s >> entries_cols;
        s >> entries_data_real;
        if (entries_flags & complex_not_real) s >> entries_data_imag;
        s >> entries_indices;
        s >> entries_indptr;
        s >> basis_flags;
        s >> basis_rows;
        s >> basis_cols;
        s >> basis_data_real;
        if (basis_flags & complex_not_real) s >> basis_data_imag;
        s >> basis_indices;
        s >> basis_indptr;

        if((((entries_flags & complex_not_real) > 0) != utils::is_complex<scalar_t>::value) ||
                (((basis_flags & complex_not_real) > 0) != utils::is_complex<scalar_t>::value)) {
            std::cout << ">>ERR" << "The data type used in the program does not fit the data type used in the serialized objects." << std::endl; // TODO throw
            abort();
        }

        // build matrix "entries_"
        std::vector<scalar_t> entries_data;
        mergeComplex(entries_data_real,entries_data_imag,entries_data);
        entries_ = eigen_sparse_t(entries_rows,entries_cols);
        entries_.makeCompressed();
        entries_.resizeNonZeros(entries_data.size());
        std::copy(entries_data.begin(),entries_data.end(),entries_.valuePtr());
        std::copy(entries_indices.begin(),entries_indices.end(),entries_.innerIndexPtr());
        std::copy(entries_indptr.begin(),entries_indptr.end(),entries_.outerIndexPtr());
        entries_.finalize();

        // build matrix "basis_"
        std::vector<scalar_t> basis_data;
        mergeComplex(basis_data_real,basis_data_imag,basis_data);
        basis_ = eigen_sparse_t(basis_rows,basis_cols);
        basis_.makeCompressed();
        basis_.resizeNonZeros(basis_data.size());
        std::copy(basis_data.begin(),basis_data.end(),basis_.valuePtr());
        std::copy(basis_indices.begin(),basis_indices.end(),basis_.innerIndexPtr());
        std::copy(basis_indptr.begin(),basis_indptr.end(),basis_.outerIndexPtr());
        basis_.finalize();
    }

    uint64_t hashEntries() {
        // TODO bring this functionality to the matrix class and use it for serialization, too
        doSerialization();
        return utils::FNV64(&bytes[0], bytes.size());
    }

    uint64_t hashBasis() {
        // TODO bring this functionality to the matrix class and use it for serialization, too
        doSerialization();
        return utils::FNV64(&bytes[0], bytes.size());
    }

    void save(std::string fname) {
        doSerialization();

        // open file
        FILE *pFile;
        pFile = fopen(fname.c_str() , "wb" );

        // write
        fwrite(&bytes[0], 1 , sizeof(byte_t)*bytes.size(), pFile );

        // close file
        fclose(pFile);
    }

    bool load(std::string fname) {
        try{
            // open file
            if (FILE *pFile = fopen(fname.c_str() , "rb" )) {
                // obtain file size:
                fseek (pFile , 0 , SEEK_END);
                size_t size_file = ftell(pFile);
                rewind(pFile);

                // read
                bytes.resize(size_file/sizeof(byte_t));
                size_t size_result = fread(&bytes[0], 1 , sizeof(byte_t)*bytes.size() , pFile);
                if (size_result != size_file) throw std::runtime_error("Matrix could not be read from file.");

                // close file
                fclose(pFile);

                doDeserialization();

                return true;
            } else {
                return false;
            }
        } catch (std::exception& e) {
#pragma omp critical(textoutput)
            std::cerr << e.what() << std::endl;
            return false;
        }
    }

    friend Hamiltonianmatrix combine(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const real_t &deltaE, std::shared_ptr<BasisnamesTwo> basis_two, const Symmetry &sym) {
        // TODO program a faster method for samebasis == true

        size_t num_basisvectors = lhs.num_basisvectors()*rhs.num_basisvectors();
        size_t num_coordinates = lhs.num_coordinates()*rhs.num_coordinates();


        ////////////////////////////////////////////////////////
        ////// Mapping used in case of reflection symmetry /////
        ////////////////////////////////////////////////////////

        std::vector<size_t> mapping(num_coordinates, -1);
        if (sym.reflection != NA) {
            std::unordered_map<StateTwo, size_t> buffer;
            for (auto state: *basis_two) {
                if (state.m[0] < 0) {
                    continue;
                }
                state.m[0] *= -1;
                state.m[1] *= -1;
                buffer[state] = state.idx;
            }
            for (auto state: *basis_two) {
                if (state.m[0] > 0) {
                    continue;
                }
                mapping[buffer[state]] = state.idx;
            }
        }


        ////////////////////////////////////////////////////////
        ////// Combine basis and entries ///////////////////////
        ////////////////////////////////////////////////////////

        // This approach does only work if lhs.entries() and rhs.entries() are diagonal // TODO assert

        // === Initialize matrices ===
        eigen_vector_t diag1 = lhs.entries().diagonal();
        eigen_vector_t diag2 = rhs.entries().diagonal();

        // Number of elements for which space sould be reserved
        size_t size_basis = num_basisvectors; // TODO estimate better
        size_t size_entries = num_basisvectors;

        Hamiltonianmatrix mat(size_basis, size_entries);

        // === Combine basis and entries ===
        size_t col = 0; // basis vector
        for (eigen_idx_t col_1=0; col_1<lhs.basis().outerSize(); ++col_1) { // outerSize() == num_cols = num_basisvectors()
            for (eigen_idx_t col_2=0; col_2<rhs.basis().outerSize(); ++col_2) {

                // In case of refelction symmetry: skip half of the basis vector pairs
                if ((sym.inversion == EVEN && col_1 <= col_2) || // gerade
                        (sym.inversion == ODD && col_1 < col_2)) { // ungerade
                    continue;
                }

                // --- Combine diagonal elements for mat.entries() ---
                scalar_t val_entries = diag1[col_1] + diag2[col_2]; // diag(V) x I + I x diag(V)

                // Check whether the new diagonal element is within the energy cutoff
                if (std::abs(val_entries) < deltaE+1e-11 || deltaE < 0) { // TODO avoid the "+1e-11" hack

                    // Variable that stores whether the combindes basis vector, that belongs to the combined diagonal element, is valid
                    bool existing = false;

                    // --- Combine basis vectors for mat.basis() ---
                    for (eigen_iterator_t triple_1(lhs.basis(),col_1); triple_1; ++triple_1) {
                        for (eigen_iterator_t triple_2(rhs.basis(),col_2); triple_2; ++triple_2) {
                            size_t row = rhs.num_coordinates()*triple_1.row() + triple_2.row(); // coordinate

                            // Get pair state that belongs to the current coordinate of the combined basis vector
                            const auto &state = basis_two->get(row);

                            float M = state.m[0]+state.m[1];
                            int parityL = std::pow(-1, state.l[0] + state.l[1]);

                            // In case of rotation symmetry: skip coordinates with wrong total magnetic momentum
                            if(sym.rotation != NA && sym.rotation != M) {
                                continue;
                            }

                            // In case of inversion and permutation symmetry: skip coordinates with wrong orbital parity
                            // It is assumed that sym.orbitalparity != NA iff inversion and permutation symmetry // TODO assert (inside Symmetry class)
                            if(sym.orbitalparity != NA && sym.orbitalparity != parityL) {
                                continue;
                            }

                            // In case of reflection symmetry: skip half of the coordinates
                            if (sym.reflection != NA && state.m[0] < 0) { // plus or minus
                                continue;
                            }

                            // Calculate coefficient that belongs to the current coordinate
                            scalar_t val_basis = triple_1.value() * triple_2.value(); // coefficient

                            // Save the coefficient taking into account the symmetrization

                            // TODO adapt code for more than one symmetry at the same time
                            // TODO handle permutation symmetry, too

                            if (sym.reflection != NA) {
                                // This code assumes the total magnetic quantum number to be zero (thus, parityM=1) // TODO assert
                                int parityJ = std::pow(-1, state.j[0] + state.j[1]);

                                val_basis /= std::sqrt(2);
                                mat.addBasis(row,col,val_basis);
                                row = mapping[row];
                                val_basis *= (sym.reflection == EVEN) ? 1*parityL*parityJ : -1*parityL*parityJ;
                            }

                            if (sym.inversion != NA && col_1 != col_2) {
                                val_basis /= std::sqrt(2);
                                mat.addBasis(row,col,val_basis);
                                row = rhs.num_coordinates()*triple_2.row() + triple_1.row();
                                val_basis *= (sym.inversion == EVEN) ? -1*parityL : 1*parityL;
                            }

                            mat.addBasis(row,col,val_basis);
                            existing = true;
                        }
                    }

                    // Save the combined diagonal element if the corresponding combined basis vector is valid
                    if (existing) {
                        mat.addEntries(col,col,val_entries);
                        ++col;
                    }
                }
            }
        }

        // Finalize the Hamiltonian matrix
        num_basisvectors = col;
        mat.compress(num_basisvectors, num_coordinates);
        return mat;
    }

    friend void energycutoff(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const real_t &deltaE, std::vector<bool> &necessary) {
        eigen_vector_t diag1 = lhs.entries().diagonal();
        eigen_vector_t diag2 = rhs.entries().diagonal();

        for (eigen_idx_t col_1=0; col_1<lhs.basis().outerSize(); ++col_1) { // outerSize() == num_cols = num_basisvectors()
            for (eigen_idx_t col_2=0; col_2<rhs.basis().outerSize(); ++col_2) {
                scalar_t val_entries = diag1[col_1] + diag2[col_2]; // diag(V) x I + I x diag(V)
                if (std::abs(val_entries) < deltaE+1e-11 || deltaE < 0) { // TODO make +1e-11 unnecessary
                    for (eigen_iterator_t triple_1(lhs.basis(),col_1); triple_1; ++triple_1) {
                        for (eigen_iterator_t triple_2(rhs.basis(),col_2); triple_2; ++triple_2) {
                            size_t row = rhs.num_coordinates()*triple_1.row() + triple_2.row(); // coordinate
                            necessary[row] = true;
                        }
                    }
                }
            }
        }
    }

protected:
    eigen_sparse_t entries_;
    eigen_sparse_t basis_;

    bytes_t bytes;

    std::vector<eigen_triplet_t> triplets_basis;
    std::vector<eigen_triplet_t> triplets_entries;
};

template <class T>
class Hamiltonian {
public:
    Hamiltonian() {}
    std::shared_ptr<Hamiltonianmatrix> get(size_t idx) {
        return matrix_diag[idx];
    }
    std::shared_ptr<const Hamiltonianmatrix> get(size_t idx) const {
        return matrix_diag[idx];
    }
    std::shared_ptr<const Configuration> getParams(size_t idx) const {
        return params[idx];
    }
    size_t size() const {
        return matrix_diag.size();
    }

    std::shared_ptr<const T> names() const {
        return basis;
    }

    void removeUnnecessaryStates(std::vector<bool> &necessary) {
        basis->removeUnnecessaryStates(necessary);
        for (auto &p : matrix_diag) {
            p->removeUnnecessaryStates(necessary);
        }
    }

protected:
    std::shared_ptr<Hamiltonianmatrix> doProcessing(std::shared_ptr<Hamiltonianmatrix> work) {
        // TODO
        return work;
    }

    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix;
    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix_diag; // TODO figure out whether the same pointers as in matrix
    std::vector<std::shared_ptr<Configuration>> params;
    std::vector<std::string> matrix_path;
    std::vector<size_t> matrix_step;
    std::vector<size_t> matrix_blocks;
    std::vector<size_t> matrix_block;
    std::vector<size_t> matrix_dimension;
    std::shared_ptr<T> basis;
};

class HamiltonianOne : public Hamiltonian<BasisnamesOne>{
public:
    HamiltonianOne(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<BasisnamesOne> basis_one) : Hamiltonian<BasisnamesOne>(), path_cache(path_cache) {
        basis = basis_one;
        configure(config);
        build();
    }



    const Configuration& getConf() const { // TODO in Configurable Klasse auslagern, von der geerbt werrden soll
        return basicconf;
    }

protected:
    void changeToSpherical(real_t val_x, real_t val_y, real_t val_z, real_t& val_p, real_t& val_m, real_t& val_0) {
        if(val_y != 0) {
            std::cout << ">>ERR" << "For fields with non-zero y-coordinates, a complex data type is needed." << std::endl;
            abort();
        }
        val_p = -val_x/std::sqrt(2);
        val_m = val_x/std::sqrt(2);
        val_0 = val_z;
    }

    void changeToSpherical(real_t val_x, real_t val_y, real_t val_z, std::complex<real_t>& val_p, std::complex<real_t>& val_m, std::complex<real_t>& val_0) {
        val_p = std::complex<real_t>(-val_x/std::sqrt(2),-val_y/std::sqrt(2));
        val_m = std::complex<real_t>(val_x/std::sqrt(2),-val_y/std::sqrt(2));
        val_0 = std::complex<real_t>(val_z,0);
    }

    void configure(const Configuration &config) {
        basicconf = basis->getConf();
        basicconf["deltaESingle"] = config["deltaESingle"];
        basicconf["diamagnetism"] = config["diamagnetism"];

        basicconf["deltaESingle"] >> deltaE;
        basicconf["species1"] >> species;

        diamagnetism = basicconf["diamagnetism"].str() == "true";

        config["minBx"] >> min_B_x;
        config["minBy"] >> min_B_y;
        config["minBz"] >> min_B_z;
        config["minEx"] >> min_E_x;
        config["minEy"] >> min_E_y;
        config["minEz"] >> min_E_z;
        config["maxBx"] >> max_B_x;
        config["maxBy"] >> max_B_y;
        config["maxBz"] >> max_B_z;
        config["maxEx"] >> max_E_x;
        config["maxEy"] >> max_E_y;
        config["maxEz"] >> max_E_z;

        if ((min_B_x == max_B_x) &&
                (min_B_y == max_B_y) &&
                (min_B_z == max_B_z) &&
                (min_E_x == max_E_x) &&
                (min_E_y == max_E_y) &&
                (min_E_z == max_E_z)) {
            nSteps = 1;
        } else {
            config["steps"] >> nSteps;
        }
    }

    void build() {
        boost::filesystem::path path_cache_mat;
        if (utils::is_complex<scalar_t>::value) {
            path_cache_mat = path_cache / "cache_matrix_complex";
        } else {
            path_cache_mat = path_cache / "cache_matrix_real";
        }
        if(!boost::filesystem::exists(path_cache_mat)){
            boost::filesystem::create_directory(path_cache_mat);
        }

        real_t tol = 1e-32;


        ////////////////////////////////////////////////////////
        ////// Build single atom basis and Hamiltonian /////////
        ////////////////////////////////////////////////////////

        // === Calculate one-atom Hamiltonian ===

        // --- Count entries of one-atom Hamiltonian ---
        size_t size_basis = basis->size();
        size_t size_energy = basis->size();

        // --- Construct one-atom  Hamiltonian and basis ---
        std::cout << "One-atom Hamiltonian, construct diagonal Hamiltonian" << std::endl;

        Hamiltonianmatrix hamiltonian_energy(size_basis, size_energy);

        real_t energy_initial = 0;
        for (const auto &state: basis->initial()) {
            energy_initial += energy_level(species,state.n,state.l,state.j);
        }
        energy_initial /= basis->initial().size(); // TODO save it to the json file

        std::vector<bool> is_necessary(basis->size(),false);
        idx_t idx = 0;
        for (const auto &state : *basis) {
            real_t val = energy_level(species,state.n,state.l,state.j)-energy_initial;
            if (std::abs(val) < deltaE+1e-11 || deltaE < 0) { // TODO
                is_necessary[state.idx] = true;
                hamiltonian_energy.addEntries(idx,idx,val);
                hamiltonian_energy.addBasis(idx,idx,1);
                ++idx;
            }
        }
        std::cout << "One-atom Hamiltonian, basis size without restrictions: " << basis->size() << std::endl;

        basis->removeUnnecessaryStates(is_necessary);

        hamiltonian_energy.compress(basis->dim(), basis->dim());

        std::cout << "One-atom Hamiltonian, basis size with restrictions: " << basis->size() << std::endl;
        std::cout << ">>BAS" << std::setw(7) << basis->size() << std::endl;

        // === Save single atom basis ===
        std::cout << "One-atom Hamiltonian, save single atom basis" << std::endl;

        // initialize uuid generator
        boost::uuids::random_generator generator;

        // generate uuid
        std::string uuid;
        boost::uuids::uuid u = generator();
        boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

        // save basis
        boost::filesystem::path path_basis = boost::filesystem::temp_directory_path();
        path_basis /= "basis_one_"+uuid+".csv";
        basis->save(path_basis.string());

        std::cout << ">>STA " << path_basis.string() << std::endl;


        ////////////////////////////////////////////////////////
        ////// Construct atom-field interaction ////////////////
        ////////////////////////////////////////////////////////

        // --- Obtain existence of fields ---
        scalar_t min_E_0, min_E_p, min_E_m, min_B_0, min_B_p, min_B_m, max_E_0, max_E_p, max_E_m, max_B_0, max_B_p, max_B_m;
        changeToSpherical(min_E_x, min_E_y, min_E_z, min_E_p, min_E_m, min_E_0);
        changeToSpherical(max_E_x, max_E_y, max_E_z, max_E_p, max_E_m, max_E_0);
        changeToSpherical(min_B_x, min_B_y, min_B_z, min_B_p, min_B_m, min_B_0);
        changeToSpherical(max_B_x, max_B_y, max_B_z, max_B_p, max_B_m, max_B_0);

        bool exist_E_0 = (std::abs(min_E_0) != 0 || std::abs(max_E_0) != 0);
        bool exist_E_1 = (std::abs(min_E_p) != 0 || std::abs(max_E_p) != 0);
        bool exist_B_0 = (std::abs(min_B_0) != 0 || std::abs(max_B_0) != 0);
        bool exist_B_1 = (std::abs(min_B_p) != 0 || std::abs(max_B_p) != 0);

        // --- Precalculate matrix elements --- // TODO parallelization
        std::cout << "One-atom Hamiltonian, precalculate matrix elements" << std::endl;

        MatrixElements matrix_elements(basicconf, species, (path_cache / "cache_elements.db").string());

        if (exist_E_0) matrix_elements.precalculateElectricMomentum(basis, 0);
        if (exist_E_1) matrix_elements.precalculateElectricMomentum(basis, 1);
        if (exist_E_1) matrix_elements.precalculateElectricMomentum(basis, -1);

        if (exist_B_0) matrix_elements.precalculateMagneticMomentum(basis, 0);
        if (exist_B_1) matrix_elements.precalculateMagneticMomentum(basis, 1);
        if (exist_B_1) matrix_elements.precalculateMagneticMomentum(basis, -1);

        if (diamagnetism && (exist_B_0 || exist_B_1)) matrix_elements.precalculateDiamagnetism(basis, 0, 0);
        if (diamagnetism && (exist_B_0 || exist_B_1)) matrix_elements.precalculateDiamagnetism(basis, 2, 0);
        if (diamagnetism && exist_B_0 && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, 1);
        if (diamagnetism && exist_B_0 && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, -1);
        if (diamagnetism && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, 2);
        if (diamagnetism && exist_B_1) matrix_elements.precalculateDiamagnetism(basis, 2, -2);

        // --- Count entries of atom-field Hamiltonian ---
        std::cout << "One-atom Hamiltonian, count number of entries within the field Hamiltonian" << std::endl;

        size_basis = basis->size();
        size_t size_electricMomentum_0 = 0;
        size_t size_electricMomentum_p = 0;
        size_t size_electricMomentum_m = 0;

        size_t size_magneticMomentum_0 = 0;
        size_t size_magneticMomentum_p = 0;
        size_t size_magneticMomentum_m = 0;

        size_t size_diamagnetism_00 = 0;
        size_t size_diamagnetism_20 = 0;
        size_t size_diamagnetism_2p = 0;
        size_t size_diamagnetism_2m = 0;
        size_t size_diamagnetism_2pp = 0;
        size_t size_diamagnetism_2mm = 0;

        for (const auto &state_col : *basis) { // TODO parallelization
            for (const auto &state_row : *basis) {
                if (state_row.idx < state_col.idx) { // lower triangle only
                    continue;
                }

                if (exist_E_0 && selectionRulesMultipole(state_row, state_col, 1, 0) ) {
                    size_electricMomentum_0++;
                } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, 1) ) {
                    size_electricMomentum_p++;
                } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, -1) ) {
                    size_electricMomentum_m++;
                }

                if (exist_B_0 &&  selectionRulesMomentum(state_row, state_col, 0) ) {
                    size_magneticMomentum_0++;
                } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, 1) ) {
                    size_magneticMomentum_p++;
                } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, -1) ) {
                    size_magneticMomentum_m++;
                }

                if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 0, 0)) {
                    size_diamagnetism_00++;
                } else if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 0)) {
                    size_diamagnetism_20++;
                } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 1)) {
                    size_diamagnetism_2p++;
                } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -1)) {
                    size_diamagnetism_2m++;
                } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 2)) {
                    size_diamagnetism_2pp++;
                } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -2)) {
                    size_diamagnetism_2mm++;
                }
            }
        }

        // --- Construct atom-field Hamiltonian ---
        std::cout << "One-atom Hamiltonian, construct field Hamiltonian" << std::endl;

        Hamiltonianmatrix hamiltonian_electricMomentum_0(size_basis, size_electricMomentum_0);
        Hamiltonianmatrix hamiltonian_electricMomentum_p(size_basis, size_electricMomentum_p);
        Hamiltonianmatrix hamiltonian_electricMomentum_m(size_basis, size_electricMomentum_m);

        Hamiltonianmatrix hamiltonian_magneticMomentum_0(size_basis, size_magneticMomentum_0);
        Hamiltonianmatrix hamiltonian_magneticMomentum_p(size_basis, size_magneticMomentum_p);
        Hamiltonianmatrix hamiltonian_magneticMomentum_m(size_basis, size_magneticMomentum_m);

        Hamiltonianmatrix hamiltonian_diamagnetism_00(size_basis, size_diamagnetism_00);
        Hamiltonianmatrix hamiltonian_diamagnetism_20(size_basis, size_diamagnetism_20);
        Hamiltonianmatrix hamiltonian_diamagnetism_2p(size_basis, size_diamagnetism_2p);
        Hamiltonianmatrix hamiltonian_diamagnetism_2m(size_basis, size_diamagnetism_2m);
        Hamiltonianmatrix hamiltonian_diamagnetism_2pp(size_basis, size_diamagnetism_2pp);
        Hamiltonianmatrix hamiltonian_diamagnetism_2mm(size_basis, size_diamagnetism_2mm);

        for (const auto &state_col : *basis) { // TODO parallelization
            for (const auto &state_row : *basis) {
                if (state_row.idx < state_col.idx) {
                    continue;
                }

                if (state_row.idx == state_col.idx) {
                    hamiltonian_electricMomentum_0.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_electricMomentum_p.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_electricMomentum_m.addBasis(state_row.idx,state_col.idx,1);

                    hamiltonian_magneticMomentum_0.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_magneticMomentum_p.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_magneticMomentum_m.addBasis(state_row.idx,state_col.idx,1);

                    hamiltonian_diamagnetism_00.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_diamagnetism_20.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_diamagnetism_2p.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_diamagnetism_2m.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_diamagnetism_2pp.addBasis(state_row.idx,state_col.idx,1);
                    hamiltonian_diamagnetism_2mm.addBasis(state_row.idx,state_col.idx,1);
                }

                if (exist_E_0 && selectionRulesMultipole(state_row, state_col, 1, 0) ) {
                    real_t val = matrix_elements.getElectricMomentum(state_row, state_col);
                    if (std::abs(val) > tol) {
                        hamiltonian_electricMomentum_0.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, 1) ) {
                    real_t val = matrix_elements.getElectricMomentum(state_row, state_col);
                    if (std::abs(val) > tol) {
                        hamiltonian_electricMomentum_p.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (exist_E_1 && selectionRulesMultipole(state_row, state_col, 1, -1) ) {
                    real_t val = matrix_elements.getElectricMomentum(state_row, state_col);
                    if (std::abs(val) > tol) {
                        hamiltonian_electricMomentum_m.addEntries(state_row.idx,state_col.idx,val);
                    }
                }

                if (exist_B_0 && selectionRulesMomentum(state_row, state_col, 0) ) {
                    real_t val = matrix_elements.getMagneticMomentum(state_row, state_col);
                    if (std::abs(val) > tol) {
                        hamiltonian_magneticMomentum_0.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, 1) ) {
                    real_t val = matrix_elements.getMagneticMomentum(state_row, state_col);
                    if (std::abs(val) > tol) {
                        hamiltonian_magneticMomentum_p.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (exist_B_1 && selectionRulesMomentum(state_row, state_col, -1) ) {
                    real_t val = matrix_elements.getMagneticMomentum(state_row, state_col);
                    if (std::abs(val) > tol) {
                        hamiltonian_magneticMomentum_m.addEntries(state_row.idx,state_col.idx,val);
                    }
                }

                if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 0, 0)) {
                    real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 0);
                    if (std::abs(val) > tol) {
                        hamiltonian_diamagnetism_00.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (diamagnetism && (exist_B_0 || exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 0)) {
                    real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                    if (std::abs(val) > tol) {
                        hamiltonian_diamagnetism_20.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 1)) {
                    real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                    if (std::abs(val) > tol) {
                        hamiltonian_diamagnetism_2p.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (diamagnetism && (exist_B_0 && exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -1)) {
                    real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                    if (std::abs(val) > tol) {
                        hamiltonian_diamagnetism_2m.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, 2)) {
                    real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                    if (std::abs(val) > tol) {
                        hamiltonian_diamagnetism_2pp.addEntries(state_row.idx,state_col.idx,val);
                    }
                } else if (diamagnetism && (exist_B_1) && selectionRulesMultipole(state_row, state_col, 2, -2)) {
                    real_t val = matrix_elements.getDiamagnetism(state_row, state_col, 2);
                    if (std::abs(val) > tol) {
                        hamiltonian_diamagnetism_2mm.addEntries(state_row.idx,state_col.idx,val);
                    }
                }
            }
        }

        std::cout << "One-atom Hamiltonian, compress field Hamiltonian" << std::endl;

        hamiltonian_electricMomentum_0.compress(basis->dim(), basis->dim());
        hamiltonian_electricMomentum_p.compress(basis->dim(), basis->dim());
        hamiltonian_electricMomentum_m.compress(basis->dim(), basis->dim());

        hamiltonian_magneticMomentum_0.compress(basis->dim(), basis->dim());
        hamiltonian_magneticMomentum_p.compress(basis->dim(), basis->dim());
        hamiltonian_magneticMomentum_m.compress(basis->dim(), basis->dim());

        hamiltonian_diamagnetism_00.compress(basis->dim(), basis->dim());
        hamiltonian_diamagnetism_20.compress(basis->dim(), basis->dim());
        hamiltonian_diamagnetism_2p.compress(basis->dim(), basis->dim());
        hamiltonian_diamagnetism_2m.compress(basis->dim(), basis->dim());
        hamiltonian_diamagnetism_2pp.compress(basis->dim(), basis->dim());
        hamiltonian_diamagnetism_2mm.compress(basis->dim(), basis->dim());


        ////////////////////////////////////////////////////////
        ////// Prepare processing of Hamiltonians //////////////
        ////////////////////////////////////////////////////////

        // TODO Put the logic in its own class

        std::cout << "One-atom Hamiltonian, processe Hamiltonians" << std::endl;

        // === Open database ===
        boost::filesystem::path path_db;

        if (utils::is_complex<scalar_t>::value) {
            path_db = path_cache / "cache_matrix_complex.db";
        } else {
            path_db = path_cache / "cache_matrix_real.db";
        }
        sqlite::handle db(path_db.string());

        // === Initialize variables ===
        bool flag_perhapsmissingtable = true;

        matrix_path.resize(nSteps);
        matrix_diag.resize(nSteps); // TODO maybe remove
        params.resize(nSteps); // TODO maybe remove


        ////////////////////////////////////////////////////////
        ////// Loop through steps //////////////////////////////
        ////////////////////////////////////////////////////////

        std::cout << ">>TOT" << std::setw(7) << nSteps << std::endl;

#pragma omp parallel for schedule(static, 1)

        // Loop through steps
        for (size_t step = 0; step < nSteps; ++step) {

            // === Get parameters for the current position inside the loop ===

            // Get fields
            real_t normalized_position = (nSteps > 1) ? step/(nSteps-1.) : 0;

            real_t Ex = min_E_x+normalized_position*(max_E_x-min_E_x);
            real_t Ey = min_E_y+normalized_position*(max_E_y-min_E_y);
            real_t Ez = min_E_z+normalized_position*(max_E_z-min_E_z);
            real_t Bx = min_B_x+normalized_position*(max_B_x-min_B_x);
            real_t By = min_B_y+normalized_position*(max_B_y-min_B_y);
            real_t Bz = min_B_z+normalized_position*(max_B_z-min_B_z);

            scalar_t E_0 = min_E_0+normalized_position*(max_E_0-min_E_0);
            scalar_t E_p = min_E_p+normalized_position*(max_E_p-min_E_p);
            scalar_t E_m = min_E_m+normalized_position*(max_E_m-min_E_m);
            scalar_t B_0 = min_B_0+normalized_position*(max_B_0-min_B_0);
            scalar_t B_p = min_B_p+normalized_position*(max_B_p-min_B_p);
            scalar_t B_m = min_B_m+normalized_position*(max_B_m-min_B_m);

            // Get configuration and save fields
            Configuration conf = basicconf;
            conf["Ex"] = Ex;
            conf["Ey"] = Ey;
            conf["Ez"] = Ez;
            conf["Bx"] = Bx;
            conf["By"] = By;
            conf["Bz"] = Bz;

            // === Create table if necessary ===
            std::stringstream query;
            std::string spacer = "";

            if (flag_perhapsmissingtable) {
                query << "CREATE TABLE IF NOT EXISTS cache_one (uuid text NOT NULL PRIMARY KEY, "
                         "created TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                         "accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP";
                for (auto p: conf) {
                    query << ", " << p.key << " text";
                }
                query << ", UNIQUE (";
                for (auto p: conf) {
                    query << spacer << p.key;
                    spacer = ", ";
                }
                query << "));";

                flag_perhapsmissingtable = false;
            }

            // === Get uuid as filename === // TODO put code in its own method
            std::string uuid = "";
            spacer = "";
            query << "SELECT uuid FROM cache_one WHERE ";
            for (auto p: conf) {
                query << spacer << p.key << "='" << p.value.str() << "'";
                spacer = " AND ";
            }
            query << ";";

#pragma omp critical(database)
            {
                sqlite::result result = db.query(query);
                if (result.size() == 1) {
                    uuid = result.first();
                }
            }

            if (uuid != "") {
                query.str(std::string());
                query << "UPDATE cache_one SET accessed = CURRENT_TIMESTAMP WHERE uuid = '" << uuid << "';";
#pragma omp critical(database)
                db.exec(query.str()); // TODO check whether this slows down the program

            } else {
                boost::uuids::uuid u = generator();
                boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

                query.str(std::string());
                query << "INSERT INTO cache_one (uuid";
                for (auto p: conf) {
                    query << ", " << p.key;
                }
                query << ") values ( '" << uuid << "'";
                for (auto p: conf) {
                    query << ", " << "'" << p.value.str() << "'";
                }
                query << ");";
#pragma omp critical(database)
                db.exec(query);
            }

            // === Check existence of files === // TODO put code in its own method

            // Check whether .mat and .json file exists and compare settings in program with settings in .json file
            boost::filesystem::path path, path_mat, path_json;

            path = path_cache_mat / ("one_" + uuid);
            path_mat = path;
            path_mat.replace_extension(".mat");
            path_json = path;
            path_json.replace_extension(".json");

            bool is_existing = false;
            if (boost::filesystem::exists(path_mat)) {
                if (boost::filesystem::exists(path_json)) {
                    Configuration params_loaded;
                    params_loaded.load_from_json(path_json.string());
                    if (conf == params_loaded) {
                        is_existing = true;
                    }
                }
            }

            // Create .json file if "is_existing" is false
            if (!is_existing) {
                conf.save_to_json(path_json.string());
            }

            // === Build and diagonalize total matrix if not existent ===
            Hamiltonianmatrix totalmatrix;

            // calculate Hamiltonian if "is_existing" is false
            std::shared_ptr<Hamiltonianmatrix> mat;
            if (!is_existing || !totalmatrix.load(path_mat.string())) {

                // --- Build matrix ---
                totalmatrix = hamiltonian_energy
                        -hamiltonian_electricMomentum_0*E_0
                        +hamiltonian_electricMomentum_p*E_m
                        +hamiltonian_electricMomentum_m*E_p
                        +hamiltonian_magneticMomentum_0*B_0
                        -hamiltonian_magneticMomentum_p*B_m
                        -hamiltonian_magneticMomentum_m*B_p
                        +hamiltonian_diamagnetism_00*(B_0*B_0-2.*B_p*B_m)
                        -hamiltonian_diamagnetism_20*(B_0*B_0+B_p*B_m)
                        +std::sqrt(3)*hamiltonian_diamagnetism_2p*B_0*B_m
                        +std::sqrt(3)*hamiltonian_diamagnetism_2m*B_0*B_p
                        -std::sqrt(1.5)*hamiltonian_diamagnetism_2pp*B_m*B_m
                        -std::sqrt(1.5)*hamiltonian_diamagnetism_2mm*B_p*B_p;

                // Stdout: Hamiltonian assembled
#pragma omp critical(textoutput)
                std::cout << ">>DIM" << std::setw(7) << totalmatrix.num_basisvectors() << std::endl
                          << "One-atom Hamiltonian, " <<  step+1 << ". Hamiltonian assembled" << std::endl;

                // --- Diagonalize matrix and save diagonalized matrix ---
                totalmatrix.diagonalize();
                totalmatrix.save(path_mat.string());

                // Stdout: Hamiltonian diagonalized
#pragma omp critical(textoutput)
                std::cout << ">>OUT" << std::setw(7) << step+1 << std::setw(7) << step << std::setw(7) << 1 << std::setw(7) << 0 << " " << path.string() << std::endl
                          << "One-atom Hamiltonian, " <<  step+1 << ". Hamiltonian diagonalized" << std::endl;
            } else {
                // Stdout: Hamiltonian loaded
#pragma omp critical(textoutput)
                std::cout << ">>DIM" << std::setw(7) << totalmatrix.num_basisvectors() << std::endl
                          << ">>OUT" << std::setw(7) << step+1 << std::setw(7) << step << std::setw(7) << 1 << std::setw(7) << 0 << " " << path.string() << std::endl
                          << "One-atom Hamiltonian, " <<  step+1 << ". Hamiltonian loaded" << std::endl;
            }

            // === Store path to configuration and diagonalized matrix ===
            matrix_path[step] = path.string();
            matrix_diag[step] = std::make_shared<Hamiltonianmatrix>(totalmatrix); // TODO maybe remove
            params[step] = std::make_shared<Configuration>(conf); // TODO maybe remove
        }

        std::cout << "One-atom Hamiltonian, all Hamiltonians processed" << std::endl;

    }

private:
    Configuration basicconf;
    real_t deltaE;
    real_t min_E_x,min_E_y,min_E_z,max_E_x,max_E_y,max_E_z,min_B_x,min_B_y,min_B_z,max_B_x,max_B_y,max_B_z;
    size_t nSteps;
    bool diamagnetism;
    std::string species;
    boost::filesystem::path path_cache;

};

class HamiltonianTwo : public Hamiltonian<BasisnamesTwo> {
public:
    HamiltonianTwo(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<HamiltonianOne> hamiltonian_one)  :
        Hamiltonian<BasisnamesTwo>(), hamiltonian_one1(hamiltonian_one), hamiltonian_one2(hamiltonian_one), path_cache(path_cache) { // TODO

        samebasis = true;

        calculate(config);
    }

    HamiltonianTwo(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<HamiltonianOne> hamiltonian_one1, std::shared_ptr<HamiltonianOne> hamiltonian_one2) :
        Hamiltonian<BasisnamesTwo>(), hamiltonian_one1(hamiltonian_one1), hamiltonian_one2(hamiltonian_one2), path_cache(path_cache) {

        samebasis = false;

        calculate(config);
    }

    void calculate(const Configuration &conf_tot) {
        boost::filesystem::path path_cache_mat;
        if (utils::is_complex<scalar_t>::value) {
            path_cache_mat = path_cache / "cache_matrix_complex";
        } else {
            path_cache_mat = path_cache / "cache_matrix_real";
        }
        if(!boost::filesystem::exists(path_cache_mat)){
            boost::filesystem::create_directory(path_cache_mat);
        }

        real_t tol = 1e-32;

        if (hamiltonian_one1->size() != hamiltonian_one2->size()) {
            std::cout << "The number of single atom Hamiltonians must be the same for both atoms." << std::endl;
            abort();
        }

        size_t nSteps_one = hamiltonian_one1->size();

        // --- generate configuration ---
        std::vector<Configuration> conf_mat;
        conf_mat.reserve(nSteps_one);

        // new, pair hamiltonian specific configuration

        if (samebasis) {
            basis = std::make_shared<BasisnamesTwo>(hamiltonian_one1->names()); // TODO remove
        } else {
            basis = std::make_shared<BasisnamesTwo>(hamiltonian_one1->names(), hamiltonian_one2->names()); // TODO remove
        }
        Configuration conf_matpair = basis->getConf();
        conf_matpair["deltaEPair"] = conf_tot["deltaEPair"];
        conf_matpair["deltaNPair"] = conf_tot["deltaNPair"];
        conf_matpair["deltaLPair"] = conf_tot["deltaLPair"];
        conf_matpair["deltaJPair"] = conf_tot["deltaJPair"];
        conf_matpair["deltaMPair"] = conf_tot["deltaMPair"];
        conf_matpair["conserveM"] = conf_tot["conserveM"];
        conf_matpair["conserveParityL"] = conf_tot["conserveParityL"];
        conf_matpair["exponent"] = conf_tot["exponent"];

        for (size_t i = 0; i < nSteps_one; ++i) {
            // old, single atom hamiltonian specific configuration
            Configuration conf_matsingle = *hamiltonian_one1->getParams(i);  // TODO
            conf_matsingle += conf_matpair;
            conf_mat.push_back(conf_matsingle);
            //conf_mat.push_back(conf_matsingle + conf_matpair); // conf_matpair overwrites settings in conf_matsingle // TODO
        }

        // setup variables
        conf_mat.back()["species1"] >> species1; // TODO order state inside cinfiguration object config.order()
        conf_mat.back()["species2"] >> species2; // TODO order state inside cinfiguration object
        conf_mat.back()["deltaEPair"] >> deltaE;
        conf_mat.back()["deltaNPair"] >> deltaN;
        conf_mat.back()["deltaLPair"] >> deltaL;
        conf_mat.back()["deltaJPair"] >> deltaJ;
        conf_mat.back()["deltaMPair"] >> deltaM;
        conserveM = conf_tot["conserveM"].str() == "true";
        conserveParityL = conf_tot["conserveParityL"].str() == "true";
        conf_tot["steps"] >> nSteps_two;
        conf_tot["minR"] >> min_R;
        conf_tot["maxR"] >> max_R;
        conf_tot["exponent"] >> multipoleexponent;

        real_t minEx, minEy, minEz, maxEx, maxEy, maxEz, minBx, minBy, minBz, maxBx, maxBy, maxBz;
        conf_tot["minEx"] >> minEx;
        conf_tot["minEy"] >> minEy;
        conf_tot["minEz"] >> minEz;
        conf_tot["maxEx"] >> maxEx;
        conf_tot["maxEy"] >> maxEy;
        conf_tot["maxEz"] >> maxEz;
        conf_tot["minBx"] >> minBx;
        conf_tot["minBy"] >> minBy;
        conf_tot["minBz"] >> minBz;
        conf_tot["maxBx"] >> maxBx;
        conf_tot["maxBy"] >> maxBy;
        conf_tot["maxBz"] >> maxBz;
        //bool fields_change_m = (minEx != 0) || (minEy != 0) || (maxEx != 0) || (maxEy != 0) || (minBx != 0) || (minBy != 0) || (maxBx != 0) || (maxBy != 0); // TODO wie richtig? so ist es eine variable, die von mehreren matrizen abhaengt
        //bool fields_change_l = (minEx != 0) || (minEy != 0) || (minEz != 0) || (maxEx != 0) || (maxEy != 0) || (maxEz != 0); // TODO wie richtig? so ist es eine variable, die von mehreren matrizen abhaengt

        //fields_change_m = true; // TODO
        //fields_change_l = true; // TODO

        if (min_R == max_R && nSteps_one == 1){
            nSteps_two = 1;
        }


        ////////////////////////////////////////////////////////
        ////// Restrict single atom states /////////////////////
        ////////////////////////////////////////////////////////

        // === Restrict states of atom 1 ===

        auto basis_one1 = hamiltonian_one1->names();
        std::vector<StateOne> initial1 = basis_one1->initial();
        std::vector<bool> necessary1(basis_one1->size(), false);

        for (const auto &state: *basis_one1) {
            bool validN = false;
            bool validL = false;
            bool validJ = false;
            bool validM = false;

            for (const auto &initial: initial1) {
                if (deltaN < 0 || std::abs(state.n - initial.n) <= deltaN) validN = true;
                if (deltaL < 0 || std::abs(state.l - initial.l) <= deltaL) validL = true;
                if (deltaJ < 0 || std::abs(state.j - initial.j) <= deltaJ) validJ = true;
                if (deltaM < 0 || std::abs(state.m - initial.m) <= deltaM) validM = true;
            }

            if (validN && validL && validJ && validM) {
                necessary1[state.idx] = true;
            }
        }

        hamiltonian_one1->removeUnnecessaryStates(necessary1);

        // === Restrict states of atom 2 ===

        if (!samebasis) {
            auto basis_one2 = hamiltonian_one2->names();
            std::vector<StateOne> initial2 = basis_one2->initial();
            std::vector<bool> necessary2(basis_one2->size(), false);

            for (const auto &state: *basis_one2) {
                bool validN = false;
                bool validL = false;
                bool validJ = false;
                bool validM = false;

                for (const auto &initial: initial2) {
                    if (deltaN < 0 || std::abs(state.n - initial.n) <= deltaN) validN = true;
                    if (deltaL < 0 || std::abs(state.l - initial.l) <= deltaL) validL = true;
                    if (deltaJ < 0 || std::abs(state.j - initial.j) <= deltaJ) validJ = true;
                    if (deltaM < 0 || std::abs(state.m - initial.m) <= deltaM) validM = true;
                }

                if (validN && validL && validJ && validM) {
                    necessary2[state.idx] = true;
                }
            }

            hamiltonian_one2->removeUnnecessaryStates(necessary2);
        }


        ////////////////////////////////////////////////////////
        ////// Build pair state basis //////////////////////////
        ////////////////////////////////////////////////////////

        // === Build pair state basis ===

        std::cout << "Two-atom Hamiltonian, build pair state basis" << std::endl;

        if (samebasis) {
            basis = std::make_shared<BasisnamesTwo>(hamiltonian_one1->names());
        } else {
            basis = std::make_shared<BasisnamesTwo>(hamiltonian_one1->names(), hamiltonian_one2->names());
        }

        std::cout << "Two-atom Hamiltonian, basis size without restrictions: " << basis->size() << std::endl;

        // === Determine necessary symmetries ===
        std::cout << "Two-atom Hamiltonian, determine symmetrized subspaces" << std::endl;
        StateTwo initial = basis->initial();

        Symmetry sym;
        sym.inversion = NA;
        sym.reflection = NA;
        sym.permutation = NA;
        sym.rotation = NA;
        sym.orbitalparity = NA;

        if (conserveM) {
            sym.rotation = initial.m[0]+initial.m[1];
        }
        if (conserveParityL) {
            sym.orbitalparity = static_cast<parity_t>(std::pow(-1, initial.l[0] + initial.l[1]));
        }

        std::vector<Symmetry> symmetries;
        if (samebasis) {
            sym.inversion = ODD;
            symmetries.push_back(sym);
            if (initial.first() != initial.second()) {
                sym.inversion = EVEN;
                symmetries.push_back(sym);
            }
        } else {
            symmetries.push_back(sym);
        }

        // TODO make use of all the symmetries

        // === Build up the list of necessary pair states ===
        std::cout << "Two-atom Hamiltonian, build up the list of necessary pair states" << std::endl;

        // Apply energy cutoff
        std::vector<bool> necessary_tmp(basis->size(), false);

#pragma omp parallel for
        for (size_t i = 0; i < nSteps_one; ++i) {
            energycutoff(*(hamiltonian_one1->get(i)), *(hamiltonian_one2->get(i)), deltaE, necessary_tmp);
        }

        // Apply restrictions due to symmetries
        std::vector<bool> necessary(basis->size(), false);

        for (const auto &state: *basis) {
            for (Symmetry sym : symmetries) {
                float M = state.m[0]+state.m[1];
                int parityL = std::pow(-1, state.l[0] + state.l[1]);

                // In case of rotation symmetry: skip pair states with wrong total magnetic momentum
                if(sym.rotation != NA && sym.rotation != M) {
                    continue;
                }

                // In case of inversion and permutation symmetry: skip pair states with wrong orbital parity
                if(sym.orbitalparity != NA && sym.orbitalparity != parityL) {
                    continue;
                }

                necessary[state.idx] = necessary_tmp[state.idx];
            }
        }

        int numNecessary = std::count(necessary.begin(), necessary.end(), true);
        std::cout << "Two-atom Hamiltonian, basis size with restrictions: " << numNecessary << std::endl;
        std::cout << ">>BAS" << std::setw(7) << numNecessary << std::endl;

        // === Save pair state basis ===
        std::cout << "Two-atom Hamiltonian, save pair state basis" << std::endl;

        // initialize uuid generator
        boost::uuids::random_generator generator;

        // generate uuid
        std::string uuid;
        boost::uuids::uuid u = generator();
        boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

        // save pair state basis
        boost::filesystem::path path_basis = boost::filesystem::temp_directory_path();
        path_basis /= "basis_two_"+uuid+".csv";
        basis->save(path_basis.string()); // TODO save only necessary entries, i.e. save pair state basis in sparse format (possibility, remove basis states but keep their idx - this would also make "if (necessary) continue" unneeded; then, "combine" has to check existence of basis element and the python script has to be adapted)

        std::cout << ">>STA " << path_basis.string() << std::endl;


        ////////////////////////////////////////////////////////
        ////// Construct atom-atom interaction /////////////////
        ////////////////////////////////////////////////////////

        // Construct pair Hamiltonians for all orders of the multipole expansion

        std::vector<int> exponent_multipole;
        std::vector<Hamiltonianmatrix> mat_multipole;
        MatrixElements matrixelements_atom1(conf_tot, species1, (path_cache / "cache_elements.db").string());
        MatrixElements matrixelements_atom2(conf_tot, species2, (path_cache / "cache_elements.db").string());
        std::vector<idx_t> size_mat_multipole;

        int idx_multipole_max = -1;

        if (multipoleexponent > 2) {

            // --- Initialize two-atom interaction Hamiltonians ---
            std::cout << "Two-atom Hamiltonian, initialize interaction Hamiltonians" << std::endl;

            int kappa_min = 1; // spherical dipole operators
            int kappa_max = multipoleexponent-kappa_min-1;
            int sumOfKappas_min = kappa_min+kappa_min;
            int sumOfKappas_max = kappa_max+kappa_min;
            idx_multipole_max = sumOfKappas_max-sumOfKappas_min;

            exponent_multipole.reserve(idx_multipole_max+1);
            mat_multipole.reserve(idx_multipole_max+1);
            size_mat_multipole.resize(idx_multipole_max+1);

            // --- Precalculate matrix elements --- // TODO parallelization
            std::cout << "Two-atom Hamiltonian, get one-atom states needed for the pair state basis"<< std::endl;

            auto basis_one1_needed = std::make_shared<BasisnamesOne>(BasisnamesOne::fromFirst(basis));
            auto basis_one2_needed = std::make_shared<BasisnamesOne>(BasisnamesOne::fromSecond(basis));

            for (int kappa = kappa_min; kappa<=kappa_max; ++kappa) {
                std::cout << "Two-atom Hamiltonian, precalculate matrix elements for kappa = " << kappa << std::endl;
                matrixelements_atom1.precalculateMultipole(basis_one1_needed, kappa);
                matrixelements_atom2.precalculateMultipole(basis_one2_needed, kappa);
            }

            // TODO if (samebasis) ...

            // --- Count entries of two-atom interaction Hamiltonians ---
            std::cout << "Two-atom Hamiltonian, count number of entries within the interaction Hamiltonians" << std::endl;

            for (int sumOfKappas = sumOfKappas_min; sumOfKappas<=sumOfKappas_max; ++sumOfKappas) {
                int idx_multipole = sumOfKappas-sumOfKappas_min;

                for (const auto &state_col : *basis) { // TODO parallelization
                    if (!necessary[state_col.idx]) continue;

                    int M_col = state_col.first().m + state_col.second().m;

                    for (const auto &state_row : *basis) {
                        if (!necessary[state_row.idx]) continue;

                        if (state_row.idx < state_col.idx) continue;
                        int M_row = state_row.first().m + state_row.second().m;
                        if (M_col != M_row) continue;

                        // multipole interaction with 1/R^(sumOfKappas+1) = 1/R^(idx_multipole+3) decay
                        for (int kappa1 = kappa_min; kappa1 <= sumOfKappas-1; ++kappa1) {
                            int kappa2 = sumOfKappas - kappa1;

                            // allowed deltaL, deltaJ, and deltaM?
                            if (selectionRulesMultipole(state_row.first(), state_col.first(), kappa1) && selectionRulesMultipole(state_row.second(), state_col.second(), kappa2)) {
                                int q1 = state_row.first().m-state_col.first().m;
                                int q2 = state_row.second().m-state_col.second().m;

                                // total momentum preserved?
                                if (q1 == -q2) {
                                    size_mat_multipole[idx_multipole]++;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // --- Construct two-atom interaction Hamiltonians ---
            size_t size_basis = basis->size();

            for (int sumOfKappas = sumOfKappas_min; sumOfKappas<=sumOfKappas_max; ++sumOfKappas) {
                std::cout << "Two-atom Hamiltonian, construct interaction Hamiltonian that belongs to 1/R^" << sumOfKappas+1 << std::endl;

                int idx_multipole = sumOfKappas-sumOfKappas_min;

                exponent_multipole.push_back(sumOfKappas+1);
                mat_multipole.push_back(Hamiltonianmatrix(size_basis, 2*size_mat_multipole[idx_multipole])); // factor of 2 because triangular matrix is not sufficient

                for (const auto &state_col : *basis) { // TODO parallelization
                    if (!necessary[state_col.idx]) continue;

                    int M_col = state_col.first().m + state_col.second().m;

                    for (const auto &state_row : *basis) {
                        if (!necessary[state_row.idx]) continue;

                        if (state_row.idx < state_col.idx) continue;
                        int M_row = state_row.first().m + state_row.second().m;
                        if (M_col != M_row) continue;

                        // construct basis
                        if (state_row.idx == state_col.idx) {
                            mat_multipole[idx_multipole].addBasis(state_row.idx,state_col.idx,1);
                        }

                        // multipole interaction with 1/R^(sumOfKappas+1) = 1/R^(idx_multipole+3) decay
                        real_t val = 0;

                        for (int kappa1 = kappa_min; kappa1 <= sumOfKappas-1; ++kappa1) {
                            int kappa2 = sumOfKappas - kappa1;

                            // allowed deltaL, deltaJ, and deltaM?
                            if (selectionRulesMultipole(state_row.first(), state_col.first(), kappa1) && selectionRulesMultipole(state_row.second(), state_col.second(), kappa2)) {
                                int q1 = state_row.first().m-state_col.first().m;
                                int q2 = state_row.second().m-state_col.second().m;

                                // total momentum preserved?
                                if (q1 == -q2) {
                                    double binomials = boost::math::binomial_coefficient<double>(kappa1+kappa2, kappa1+q1)*boost::math::binomial_coefficient<double>(kappa1+kappa2, kappa2-q2);
                                    val += std::pow(-1,kappa2) * std::sqrt(binomials) * matrixelements_atom1.getMultipole(state_row.first(), state_col.first(), kappa1)*
                                            matrixelements_atom2.getMultipole(state_row.second(), state_col.second(), kappa2);
                                }
                            }
                        }

                        if (std::abs(val) > tol) {
                            mat_multipole[idx_multipole].addEntries(state_row.idx,state_col.idx,val);
                            if (state_row.idx != state_col.idx) mat_multipole[idx_multipole].addEntries(state_col.idx,state_row.idx,val); // triangular matrix is not sufficient because of basis change
                        }

                        // TODO state_two soll std::array<state_one, 2> sein! Dann geht auch die Abfrage der selection rules eindeutiger
                    }
                }

                std::cout << "Two-atom Hamiltonian, compress interaction Hamiltonian that belongs to 1/R^" << sumOfKappas+1 << std::endl;

                mat_multipole[idx_multipole].compress(basis->dim(), basis->dim()); // TODO substitute dim() by size()
            }
        }


        ////////////////////////////////////////////////////////
        ////// Prepare processing of Hamiltonians //////////////
        ////////////////////////////////////////////////////////

        // TODO Put the logic in its own class

        std::cout << "Two-atom Hamiltonian, process Hamiltonians" << std::endl;

        // === Open database ===
        boost::filesystem::path path_db;

        if (utils::is_complex<scalar_t>::value) {
            path_db = path_cache / "cache_matrix_complex.db";
        } else {
            path_db = path_cache / "cache_matrix_real.db";
        }
        sqlite::handle db(path_db.string());

        // === Initialize variables ===
        bool flag_perhapsmissingtable = true;

        std::map<parity_t,std::string> symmetries_name;
        symmetries_name[EVEN] = "sym";
        symmetries_name[ODD] = "asym";
        symmetries_name[NA] = "all";

        matrix_path.resize(nSteps_two*symmetries.size());

        // --- Determine combined single atom matrices ---
        // Construct pair Hamiltonian consistent of combined one-atom Hamiltonians (1 x Hamiltonian2 + Hamiltonian1 x 1)

        std::vector<Hamiltonianmatrix> mat_single;

        // Check if one_atom Hamiltonians change with step_two
        // It is assumed that nSteps_one = 1 if nSteps_two != nSteps_one // TODO introduce variable "is_mat_single_const" to improve readability
        if (nSteps_two != nSteps_one) {
            std::cout << "Two-atom Hamiltonian, construct contribution of combined one-atom Hamiltonians" << std::endl;

            mat_single.resize(symmetries.size());

#pragma omp parallel for
            for (size_t idx_symmetry = 0; idx_symmetry < symmetries.size(); ++idx_symmetry) {
                Symmetry sym = symmetries[idx_symmetry];

                // Combine the Hamiltonians of the two atoms
                mat_single[idx_symmetry] = combine(*(hamiltonian_one1->get(0)), *(hamiltonian_one2->get(0)), deltaE, basis, sym);

                // Remove more or less empty basis vectors
                mat_single[idx_symmetry].removeUnnecessaryBasisvectors();
            }
        }

        // --- Determine transformed interaction matrices ---
        std::vector<Hamiltonianmatrix> mat_multipole_transformed;

        // Check if one_atom Hamiltonians change with step_two
        if (nSteps_two != nSteps_one) {
            std::cout << "Two-atom Hamiltonian, construct transformed interaction matrices" << std::endl;

            mat_multipole_transformed.resize(symmetries.size()*(idx_multipole_max+1));

#pragma omp parallel for collapse(2)
            for (size_t idx_symmetry = 0; idx_symmetry < symmetries.size(); ++idx_symmetry) {
                for (int idx_multipole = 0; idx_multipole <= idx_multipole_max; ++idx_multipole) {
                    mat_multipole_transformed[idx_symmetry*(idx_multipole_max+1)+idx_multipole] = mat_multipole[idx_multipole].changeBasis(mat_single[idx_symmetry].basis());
                }
            }
        }


        ////////////////////////////////////////////////////////
        ////// Loop through steps and symmetries ///////////////
        ////////////////////////////////////////////////////////

        std::cout << ">>TOT" << std::setw(7) << nSteps_two*symmetries.size() << std::endl;

#pragma omp parallel for collapse(2) schedule(static, 1)

        // Loop through steps
        for (size_t step_two = 0; step_two < nSteps_two; ++step_two) {

            // Loop through symmetries
            for (size_t idx_symmetry = 0; idx_symmetry < symmetries.size(); ++idx_symmetry) {
                Symmetry sym = symmetries[idx_symmetry];

                size_t step = step_two*symmetries.size()+idx_symmetry;

                // === Get parameters for the current position inside the loop ===
                int single_idx = (nSteps_two == nSteps_one) ? step_two : 0;

                // Get interatomic distance
                real_t normalized_position = (nSteps_two > 1) ? step_two/(nSteps_two-1.) : 0;
                real_t position = min_R+normalized_position*(max_R-min_R);

                // Get configuration and save postions and symmetries
                Configuration conf = conf_mat[single_idx];
                conf["R"] = position;
                conf["symmetry"] = symmetries_name[sym.inversion]; // TODO adapt for other symmetries
                conf["sub"] = 0; // TODO remove

                // === Create table if necessary ===
                std::stringstream query;
                std::string spacer = "";

                if (flag_perhapsmissingtable) {
                    query << "CREATE TABLE IF NOT EXISTS cache_two (uuid text NOT NULL PRIMARY KEY, "
                             "created TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                             "accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP";
                    for (auto p: conf) {
                        query << ", " << p.key << " text";
                    }
                    query << ", UNIQUE (";
                    for (auto p: conf) {
                        query << spacer << p.key;
                        spacer = ", ";
                    }
                    query << "));";

                    flag_perhapsmissingtable = false;
                }

                // === Get uuid as filename ===
                std::string uuid = "";
                spacer = "";
                query << "SELECT uuid FROM cache_two WHERE ";
                for (auto p: conf) {
                    query << spacer << p.key << "='" << p.value.str() << "'";
                    spacer = " AND ";
                }
                query << ";";

#pragma omp critical(database)
                {
                    sqlite::result result = db.query(query);
                    if (result.size() == 1) {
                        uuid = result.first();
                    }
                }

                if (uuid != "") {
                    query.str(std::string());
                    query << "UPDATE cache_two SET accessed = CURRENT_TIMESTAMP WHERE uuid = '" << uuid << "';";
#pragma omp critical(database)
                    db.exec(query.str()); // TODO check whether this slows down the program

                } else {
                    boost::uuids::uuid u = generator();
                    boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

                    query.str(std::string());
                    query << "INSERT INTO cache_two (uuid";
                    for (auto p: conf) {
                        query << ", " << p.key;
                    }
                    query << ") values ( '" << uuid << "'";
                    for (auto p: conf) {
                        query << ", " << "'" << p.value.str() << "'";
                    }
                    query << ");";
#pragma omp critical(database)
                    db.exec(query);
                }

                // === Check existence of files ===

                // Check whether .mat and .json file exists and compare settings in program with settings in .json file
                boost::filesystem::path path, path_mat, path_json;

                path = path_cache_mat / ("two_" + uuid);
                path_mat = path;
                path_mat.replace_extension(".mat");
                path_json = path;
                path_json.replace_extension(".json");

                bool is_existing = false;
                if (boost::filesystem::exists(path_mat)) {
                    if (boost::filesystem::exists(path_json)) {
                        Configuration params_loaded;
                        params_loaded.load_from_json(path_json.string());
                        if (conf == params_loaded) {
                            is_existing = true;
                        }
                    }
                }

                // Create .json file if "is_existing" is false
                if (!is_existing) {
                    conf.save_to_json(path_json.string());
                }

                // === Build and diagonalize total matrix if not existent ===
                Hamiltonianmatrix totalmatrix;

                if (!is_existing || !totalmatrix.load(path_mat.string())) {

                    // --- Combine single atom matrices ---
                    if (nSteps_two == nSteps_one) {
                        totalmatrix = combine(*(hamiltonian_one1->get(step_two)), *(hamiltonian_one2->get(step_two)), deltaE, basis, sym);
                        totalmatrix.removeUnnecessaryBasisvectors();
                    } else {
                        totalmatrix = mat_single[idx_symmetry];
                    }

                    // --- Add interaction ---
                    for (int idx_multipole = 0; idx_multipole <= idx_multipole_max; ++idx_multipole) {
                        real_t pos = 1./std::pow(position,exponent_multipole[idx_multipole]);
                        if (nSteps_two == nSteps_one) {
                            totalmatrix += mat_multipole[idx_multipole].changeBasis(totalmatrix.basis())*pos;
                        } else {
                            totalmatrix += mat_multipole_transformed[idx_symmetry*(idx_multipole_max+1)+idx_multipole]*pos;
                        }
                    }

                    // Stdout: Hamiltonian assembled
#pragma omp critical(textoutput)
                    std::cout << ">>DIM" << std::setw(7) << totalmatrix.num_basisvectors() << std::endl
                              << "Two-atom Hamiltonian, " <<  step+1 << ". Hamiltonian assembled" << std::endl;

                    // --- Diagonalize matrix and save diagonalized matrix ---
                    totalmatrix.diagonalize();
                    totalmatrix.save(path_mat.string());

                    // Stdout: Hamiltonian diagonalized
#pragma omp critical(textoutput)
                    std::cout << ">>OUT" << std::setw(7) << step+1 << std::setw(7) << step_two << std::setw(7) << symmetries.size() << std::setw(7) << idx_symmetry << " " << path.string() << std::endl
                              << "Two-atom Hamiltonian, " <<  step+1 << ". Hamiltonian diagonalized" << std::endl;
                } else {

                    // Stdout: Hamiltonian loaded
#pragma omp critical(textoutput)
                    std::cout << ">>DIM" << std::setw(7) << totalmatrix.num_basisvectors() << std::endl
                              << ">>OUT" << std::setw(7) << step+1 << std::setw(7) << step_two << std::setw(7) << symmetries.size() << std::setw(7) << idx_symmetry << " " << path.string() << std::endl
                              << "Two-atom Hamiltonian, " <<  step+1 << ". Hamiltonian loaded" << std::endl;
                }

                // === Store path to configuration and diagonalized matrix ===
                matrix_path[step] = path.string();
            }
        }

        std::cout << "Two-atom Hamiltonian, all Hamiltonians processed" << std::endl;
    }

private:
    std::shared_ptr<HamiltonianOne> hamiltonian_one1; // TODO const HamiltonianOne
    std::shared_ptr<HamiltonianOne> hamiltonian_one2;
    real_t deltaE;
    int deltaN;
    int deltaL;
    int deltaJ;
    int deltaM;
    size_t nSteps_two;
    std::string species1, species2;
    real_t min_R, max_R;
    int multipoleexponent;
    bool samebasis;
    bool conserveM;
    bool conserveParityL;
    boost::filesystem::path path_cache;
};


///############################################################################
///### MAIN LOOP ##############################################################
///############################################################################

int main(int argc, char **argv) {
    std::cout << std::unitbuf;

    Eigen::setNbThreads(1); // TODO set it to setNbThreads(0) when Eigen's multithreading is needed

    // === Parse command line ===
    namespace po = boost::program_options;

    po::options_description desc("Usage");
    desc.add_options()
            ("help,?", "produce this help message")
            ("config,c", po::value<std::string>()->required(),"Path to config JSON file")
            ("output,o", po::value<std::string>()->required(),"Path to cache JSON file")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if ( vm.count("help") )
    {
        std::cout << desc << std::endl;
        return 0;
    }

    try
    {
        po::notify(vm);
    }
    catch (po::required_option& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    boost::filesystem::path path_config = boost::filesystem::absolute(vm["config"].as<std::string>());
    boost::filesystem::path path_cache  = boost::filesystem::absolute(vm["output"].as<std::string>());

    // === Load configuration ===
    Configuration config;
    config.load_from_json(path_config.string());

    bool existAtom1 = config.count("species1") && config.count("n1") && config.count("l1") && config.count("j1") && config.count("m1");
    bool existAtom2 = config.count("species2") && config.count("n2") && config.count("l2") && config.count("j2") && config.count("m2");

    // === Solve the system ===
    bool combined = config["samebasis"].str() == "true";

    if (combined) {
        if (config["species1"].str() != config["species2"].str()) {
            std::cout << "species1 and species2 has to be the same in order to use the same basis set." << std::endl;
            return 1;
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one;
        if (existAtom1 && existAtom2) {
            std::cout << ">>TYP" << std::setw(7) << 3 << std::endl;
            auto basisnames_one = std::make_shared<BasisnamesOne>(BasisnamesOne::fromBoth(config));
            hamiltonian_one = std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && config.count("minR")) {
            std::cout << ">>TYP" << std::setw(7) << 2 << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, hamiltonian_one);
        }
    } else {
        std::shared_ptr<HamiltonianOne> hamiltonian_one1;
        if (existAtom1) {
            std::cout << ">>TYP" << std::setw(7) << 0 << std::endl;
            auto basisnames_one1 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromFirst(config));
            hamiltonian_one1 = std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one1);
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one2;
        if (existAtom2) {
            std::cout << ">>TYP" << std::setw(7) << 1 << std::endl;
            auto basisnames_one2 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromSecond(config));
            hamiltonian_one2 = std::make_shared<HamiltonianOne>(config, path_cache, basisnames_one2);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && config.count("minR")) {
            std::cout << ">>TYP" << std::setw(7) << 2 << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, hamiltonian_one1, hamiltonian_one2);
        }
    }

    // === Communicate that everything has finished ===
    std::cout << ">>END" << std::endl;

    return 0;
}
