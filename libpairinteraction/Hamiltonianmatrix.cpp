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

#include "Hamiltonianmatrix.h"
#include <stdexcept>

#include <boost/format.hpp>

Hamiltonianmatrix::Hamiltonianmatrix() = default;

Hamiltonianmatrix::Hamiltonianmatrix(const eigen_sparse_t &entries, const eigen_sparse_t &basis)
    : entries_(entries), basis_(basis) {}

Hamiltonianmatrix::Hamiltonianmatrix(size_t szBasis, size_t szEntries) {
    triplets_basis.reserve(szBasis);
    triplets_entries.reserve(szEntries);
}

eigen_sparse_t &Hamiltonianmatrix::entries() {
    bytes.clear();
    return entries_;
}

const eigen_sparse_t &Hamiltonianmatrix::entries() const { return entries_; }

eigen_sparse_t &Hamiltonianmatrix::basis() {
    bytes.clear();
    return basis_;
}

const eigen_sparse_t &Hamiltonianmatrix::basis() const { return basis_; }

size_t Hamiltonianmatrix::num_basisvectors() const { return basis_.cols(); }

size_t Hamiltonianmatrix::num_coordinates() const { return basis_.rows(); }

void Hamiltonianmatrix::addBasis(idx_t row, idx_t col, scalar_t val) {
    triplets_basis.emplace_back(row, col, val);
}

void Hamiltonianmatrix::addEntries(idx_t row, idx_t col, scalar_t val) {
    triplets_entries.emplace_back(row, col, val);
}

void Hamiltonianmatrix::compress(size_t nBasis, size_t nCoordinates) {
    basis_.resize(nCoordinates, nBasis);
    entries_.resize(nBasis, nBasis);
    basis_.setFromTriplets(triplets_basis.begin(), triplets_basis.end());
    entries_.setFromTriplets(triplets_entries.begin(), triplets_entries.end());
    triplets_basis.clear();
    triplets_entries.clear();
}

std::vector<Hamiltonianmatrix> Hamiltonianmatrix::findSubs() const { // TODO
    std::vector<Hamiltonianmatrix> submatrices;
    submatrices.push_back(*this);
    return submatrices;
}

Hamiltonianmatrix Hamiltonianmatrix::abs() const {
    return Hamiltonianmatrix(entries_.cwiseAbs().cast<scalar_t>(), basis_);
}

Hamiltonianmatrix Hamiltonianmatrix::changeBasis(const eigen_sparse_t &basis) const {
    auto transformator = basis_.adjoint() * basis;
    auto entries = transformator.adjoint() * entries_ * transformator;
    return Hamiltonianmatrix(entries, basis);
}

void Hamiltonianmatrix::applyCutoff(double cutoff) {
    bytes.clear();

    // build transformator
    eigen_vector_t diag = entries_.diagonal();

    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(num_basisvectors());

    size_t idxBasis = 0;
    for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
        if (std::abs(diag[idx]) < cutoff) {
            triplets_transformator.emplace_back(idx, idxBasis++, 1);
        }
    }

    eigen_sparse_t transformator(this->num_basisvectors(), idxBasis);
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = basis_ * transformator;
    entries_ = transformator.adjoint() * entries_ * transformator;
}

void Hamiltonianmatrix::findUnnecessaryStates(std::vector<bool> &isNecessaryCoordinate) const {
    std::vector<double> isNecessaryCoordinate_real(num_coordinates(), 0);
    for (int k = 0; k < basis_.outerSize(); ++k) {
        for (eigen_iterator_t triple(basis_, k); triple; ++triple) {
            isNecessaryCoordinate_real[triple.row()] += std::pow(std::abs(triple.value()), 2);
        }
    }

    for (size_t idx = 0; idx < this->num_coordinates(); ++idx) {
        if (isNecessaryCoordinate_real[idx] > 0.05) { // TODO
            isNecessaryCoordinate[idx] = true;
        }
    }
}

void Hamiltonianmatrix::removeUnnecessaryBasisvectors(
    const std::vector<bool> &isNecessaryCoordinate) {
    bytes.clear();

    // build transformator
    std::vector<double> isNecessaryBasisvector(num_basisvectors(), 0);
    for (int k_1 = 0; k_1 < basis_.outerSize(); ++k_1) {
        for (eigen_iterator_t triple(basis_, k_1); triple; ++triple) {
            ptrdiff_t col = triple.col(); // basis vector
            ptrdiff_t row = triple.row(); // coordinate
            if (isNecessaryCoordinate[row]) {
                isNecessaryBasisvector[col] += std::pow(std::abs(triple.value()), 2);
            }
        }
    }

    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(num_basisvectors());

    size_t idxBasis = 0;
    for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
        if (isNecessaryBasisvector[idx] > 0.05) { // TODO
            triplets_transformator.emplace_back(idx, idxBasis++, 1);
        }
    }

    eigen_sparse_t transformator(this->num_basisvectors(), idxBasis);
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = basis_ * transformator;
    entries_ = transformator.adjoint() * entries_ * transformator;
}

void Hamiltonianmatrix::removeUnnecessaryBasisvectors() {
    bytes.clear();

    // build transformator
    std::vector<double> isNecessaryBasisvector(num_basisvectors(), 0);
    for (int k_1 = 0; k_1 < basis_.outerSize(); ++k_1) {
        for (eigen_iterator_t triple(basis_, k_1); triple; ++triple) {
            ptrdiff_t col = triple.col(); // basis vector
            isNecessaryBasisvector[col] += std::pow(std::abs(triple.value()), 2);
        }
    }

    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(num_basisvectors());

    size_t idxBasis = 0;
    for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
        if (isNecessaryBasisvector[idx] > 0.05) {
            triplets_transformator.emplace_back(idx, idxBasis++, 1);
        }
    }

    eigen_sparse_t transformator(this->num_basisvectors(), idxBasis);
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = basis_ * transformator;
    entries_ = transformator.adjoint() * entries_ * transformator;
}

void Hamiltonianmatrix::removeUnnecessaryStates(const std::vector<bool> &isNecessaryCoordinate) {
    bytes.clear();

    // build transformator
    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(num_coordinates());

    size_t idxCoordinate = 0;
    for (size_t idx = 0; idx < this->num_coordinates(); ++idx) {
        if (isNecessaryCoordinate[idx]) {
            triplets_transformator.emplace_back(idxCoordinate++, idx, 1);
        }
    }

    eigen_sparse_t transformator(idxCoordinate, this->num_coordinates());
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = transformator * basis_;
}

Hamiltonianmatrix Hamiltonianmatrix::getBlock(const std::vector<ptrdiff_t> &indices) {
    std::vector<eigen_triplet_t> triplets_transformator;
    triplets_transformator.reserve(indices.size());
    for (size_t idx = 0; idx < indices.size(); ++idx) {
        triplets_transformator.emplace_back(indices[idx], idx, 1);
    }
    eigen_sparse_t transformator(this->num_basisvectors(), indices.size());
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    eigen_sparse_t block_entries = transformator.adjoint() * entries_ * transformator;
    eigen_sparse_t block_basis = basis_ * transformator;

    return Hamiltonianmatrix(block_entries, block_basis);
}

void Hamiltonianmatrix::diagonalize() {
    if (this->num_basisvectors() > 1) { // NOLINT
        // diagonalization
        Eigen::SelfAdjointEigenSolver<eigen_dense_t> eigensolver(eigen_dense_t(this->entries()));

        // eigenvalues and eigenvectors
        eigen_vector_double_t evals = eigensolver.eigenvalues();
        eigen_sparse_t evecs = eigensolver.eigenvectors().sparseView(1e-4, 0.5);

        this->entries().setZero();
        this->entries().reserve(evals.size());
        for (int idx = 0; idx < evals.size(); ++idx) {
            this->entries().insert(idx, idx) = evals.coeffRef(idx);
        }
        this->entries().makeCompressed();

        this->basis() = (this->basis() * evecs).pruned(1e-4, 0.5);
    }
}

Hamiltonianmatrix operator+(Hamiltonianmatrix lhs, const Hamiltonianmatrix &rhs) {
    lhs.bytes.clear();
    lhs.entries_ += rhs.entries_;
    return lhs;
}

Hamiltonianmatrix operator-(Hamiltonianmatrix lhs, const Hamiltonianmatrix &rhs) {
    lhs.bytes.clear();
    lhs.entries_ -= rhs.entries_;
    return lhs;
}

Hamiltonianmatrix operator*(const scalar_t &lhs, Hamiltonianmatrix rhs) {
    rhs.bytes.clear();
    rhs.entries_ *= lhs;
    return rhs;
}

Hamiltonianmatrix operator*(Hamiltonianmatrix lhs, const scalar_t &rhs) {
    lhs.bytes.clear();
    lhs.entries_ *= rhs;
    return lhs;
}

Hamiltonianmatrix &Hamiltonianmatrix::operator+=(const Hamiltonianmatrix &rhs) {
    bytes.clear();
    entries_ += rhs.entries_;
    return *this;
}

Hamiltonianmatrix &Hamiltonianmatrix::operator-=(const Hamiltonianmatrix &rhs) {
    bytes.clear();
    entries_ -= rhs.entries_;
    return *this;
}

bytes_t &Hamiltonianmatrix::serialize() {
    doSerialization();
    return bytes;
}

void Hamiltonianmatrix::doSerialization() {
    if (bytes.empty()) {
        entries_.makeCompressed();
        basis_.makeCompressed();

        // convert matrix "entries" to vectors of primitive data types
        byte_t entries_flags = 0;
        if (entries_.IsRowMajor != 0) {
            entries_flags |= csr_not_csc;
        }
        if (utils::is_complex<scalar_t>::value) {
            entries_flags |= complex_not_real;
        }
        storage_idx_t entries_rows = entries_.rows();
        storage_idx_t entries_cols = entries_.cols();
        std::vector<scalar_t> entries_data(entries_.valuePtr(),
                                           entries_.valuePtr() + entries_.nonZeros());
        std::vector<storage_double> entries_data_real, entries_data_imag;
        splitComplex(entries_data_real, entries_data_imag, entries_data);
        std::vector<storage_idx_t> entries_indices(entries_.innerIndexPtr(),
                                                   entries_.innerIndexPtr() + entries_.nonZeros());
        std::vector<storage_idx_t> entries_indptr(entries_.outerIndexPtr(),
                                                  entries_.outerIndexPtr() + entries_.outerSize());

        // convert matrix "basis" to vectors of primitive data types
        byte_t basis_flags = 0;
        if (basis_.IsRowMajor != 0) {
            basis_flags |= csr_not_csc;
        }
        if (utils::is_complex<scalar_t>::value) {
            basis_flags |= complex_not_real;
        }
        storage_idx_t basis_rows = basis_.rows();
        storage_idx_t basis_cols = basis_.cols();
        std::vector<scalar_t> basis_data(basis_.valuePtr(), basis_.valuePtr() + basis_.nonZeros());
        std::vector<storage_double> basis_data_real, basis_data_imag;
        splitComplex(basis_data_real, basis_data_imag, basis_data);
        std::vector<storage_idx_t> basis_indices(basis_.innerIndexPtr(),
                                                 basis_.innerIndexPtr() + basis_.nonZeros());
        std::vector<storage_idx_t> basis_indptr(basis_.outerIndexPtr(),
                                                basis_.outerIndexPtr() + basis_.outerSize());

        // serialize vectors of primitive data types
        Serializer s;
        s << entries_flags;
        s << entries_rows;
        s << entries_cols;
        s << entries_data_real;
        if ((entries_flags & complex_not_real) != 0) {
            s << entries_data_imag;
        }
        s << entries_indices;
        s << entries_indptr;
        s << basis_flags;
        s << basis_rows;
        s << basis_cols;
        s << basis_data_real;
        if ((basis_flags & complex_not_real) != 0) {
            s << basis_data_imag;
        }
        s << basis_indices;
        s << basis_indptr;
        s.save(bytes);
    }
}

void Hamiltonianmatrix::deserialize(bytes_t &bytesin) {
    bytes = bytesin;
    doDeserialization();
}

void Hamiltonianmatrix::doDeserialization() {
    // deserialize vectors of primitive data types
    byte_t entries_flags;
    storage_idx_t entries_rows, entries_cols;
    std::vector<storage_double> entries_data_real, entries_data_imag;
    std::vector<idx_t> entries_indices;
    std::vector<idx_t> entries_indptr;
    byte_t basis_flags;
    storage_idx_t basis_rows, basis_cols;
    std::vector<storage_double> basis_data_real, basis_data_imag;
    std::vector<idx_t> basis_indices;
    std::vector<idx_t> basis_indptr;

    Serializer s;
    s.load(bytes);
    s >> entries_flags;
    s >> entries_rows;
    s >> entries_cols;
    s >> entries_data_real;
    if ((entries_flags & complex_not_real) != 0) {
        s >> entries_data_imag;
    }
    s >> entries_indices;
    s >> entries_indptr;
    s >> basis_flags;
    s >> basis_rows;
    s >> basis_cols;
    s >> basis_data_real;
    if ((basis_flags & complex_not_real) != 0) {
        s >> basis_data_imag;
    }
    s >> basis_indices;
    s >> basis_indptr;

    if ((((entries_flags & complex_not_real) > 0) != utils::is_complex<scalar_t>::value) ||
        (((basis_flags & complex_not_real) > 0) != utils::is_complex<scalar_t>::value)) {
        std::string msg("The data type used in the program does not fit the data type used in the "
                        "serialized objects.");
        std::cout << boost::format(">>ERR%s") % msg.c_str() << std::endl;
        throw std::runtime_error(msg);
    }

    // build matrix "entries_"
    std::vector<scalar_t> entries_data;
    mergeComplex(entries_data_real, entries_data_imag, entries_data);
    entries_ = eigen_sparse_t(entries_rows, entries_cols);
    entries_.makeCompressed();
    entries_.resizeNonZeros(entries_data.size());
    std::copy(entries_data.begin(), entries_data.end(), entries_.valuePtr());
    std::copy(entries_indices.begin(), entries_indices.end(), entries_.innerIndexPtr());
    std::copy(entries_indptr.begin(), entries_indptr.end(), entries_.outerIndexPtr());
    entries_.finalize();

    // build matrix "basis_"
    std::vector<scalar_t> basis_data;
    mergeComplex(basis_data_real, basis_data_imag, basis_data);
    basis_ = eigen_sparse_t(basis_rows, basis_cols);
    basis_.makeCompressed();
    basis_.resizeNonZeros(basis_data.size());
    std::copy(basis_data.begin(), basis_data.end(), basis_.valuePtr());
    std::copy(basis_indices.begin(), basis_indices.end(), basis_.innerIndexPtr());
    std::copy(basis_indptr.begin(), basis_indptr.end(), basis_.outerIndexPtr());
    basis_.finalize();
}

uint64_t Hamiltonianmatrix::hashEntries() {
    // TODO bring this functionality to the matrix class and use it for serialization, too
    doSerialization();
    return utils::FNV64(&bytes[0], bytes.size());
}

uint64_t Hamiltonianmatrix::hashBasis() {
    // TODO bring this functionality to the matrix class and use it for serialization, too
    doSerialization();
    return utils::FNV64(&bytes[0], bytes.size());
}

void Hamiltonianmatrix::save(const std::string &fname) {
    doSerialization();

    // open file
    FILE *pFile;
    pFile = fopen(fname.c_str(), "wb");

    // write
    fwrite(&bytes[0], 1, sizeof(byte_t) * bytes.size(), pFile);

    // close file
    fclose(pFile);
}

bool Hamiltonianmatrix::load(const std::string &fname) {
    try {
        // open file
        if (FILE *pFile = fopen(fname.c_str(), "rb")) {
            // obtain file size:
            fseek(pFile, 0, SEEK_END);
            size_t size_file = ftell(pFile);
            rewind(pFile);

            // read
            bytes.resize(size_file / sizeof(byte_t));
            size_t size_result = fread(&bytes[0], 1, sizeof(byte_t) * bytes.size(), pFile);
            if (size_result != size_file) {
                throw std::runtime_error("Matrix could not be read from file.");
            }

            // close file
            fclose(pFile);

            doDeserialization();

            return true;
        }
        return false;

    } catch (std::exception &e) {
#pragma omp critical(textoutput)
        std::cerr << e.what() << std::endl;
        return false;
    }
}

Hamiltonianmatrix combine(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs,
                          const double &deltaE, const std::shared_ptr<BasisnamesTwo> &basis_two,
                          const Symmetry &sym) {
    // TODO program a faster method for samebasis == true

    size_t num_basisvectors = lhs.num_basisvectors() * rhs.num_basisvectors();
    size_t num_coordinates = lhs.num_coordinates() * rhs.num_coordinates();

    ////////////////////////////////////////////////////////
    ////// Mapping used in case of reflection symmetry /////
    ////////////////////////////////////////////////////////

    std::vector<size_t> mapping(num_coordinates, -1);
    if (sym.reflection != NA) { // NOLINT
        std::unordered_map<StateTwo, size_t> buffer;
        for (auto state : *basis_two) {
            if (state.m[0] < 0) {
                continue;
            }
            state.m[0] *= -1;
            state.m[1] *= -1;
            buffer[state] = state.idx;
        }
        for (auto state : *basis_two) {
            if (state.m[0] > 0) {
                continue;
            }
            mapping[buffer[state]] = state.idx;
            if (sym.inversion != NA || sym.permutation != NA) {
                mapping[state.idx] = buffer[state];
            }
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
    for (int col_1 = 0; col_1 < lhs.basis().outerSize();
         ++col_1) { // outerSize() == num_cols = num_basisvectors()
        for (int col_2 = 0; col_2 < rhs.basis().outerSize(); ++col_2) {

            // In case of inversion symmetry: skip half of the basis vector pairs
            if ((sym.inversion == EVEN && col_1 <= col_2) || // gerade
                (sym.inversion == ODD && col_1 < col_2)) {   // ungerade
                continue;
            }

            // In case of permutation symmetry: skip half of the basis vector pairs
            if ((sym.permutation == EVEN && col_1 <= col_2) || // sym
                (sym.permutation == ODD && col_1 < col_2)) {   // asym
                continue;
            }

            // TODO combine inversion and permutation symmetry (one variable "permuinversion"),
            // think about further simplifications

            // --- Combine diagonal elements for mat.entries() ---
            scalar_t val_entries = diag1[col_1] + diag2[col_2]; // diag(V) x I + I x diag(V)

            // Check whether the new diagonal element is within the energy cutoff
            if (std::abs(val_entries) < deltaE + 1e-11 ||
                deltaE < 0) { // TODO avoid the "+1e-11" hack

                // Variable that stores whether the combindes basis vector, that belongs to the
                // combined diagonal element, is valid
                bool existing = false;

                // --- Combine basis vectors for mat.basis() ---
                for (eigen_iterator_t triple_1(lhs.basis(), col_1); triple_1; ++triple_1) {
                    for (eigen_iterator_t triple_2(rhs.basis(), col_2); triple_2; ++triple_2) {
                        size_t row =
                            rhs.num_coordinates() * triple_1.row() + triple_2.row(); // coordinate

                        // Get pair state that belongs to the current coordinate of the combined
                        // basis vector
                        const auto &state = basis_two->get(row);

                        float M = state.m[0] + state.m[1];
                        int parityL = std::pow(-1, state.l[0] + state.l[1]);
                        int parityJ = std::pow(-1, state.j[0] + state.j[1]);
                        int parityM = std::pow(-1, state.m[0] + state.m[1]);

                        // In case of inversion and reflection symmetry: check whether the inversion
                        // symmetric state is already reflection symmetric
                        bool skip_reflection = false;
                        if (sym.inversion != NA && col_1 != col_2 && sym.reflection != NA &&
                            mapping[row] ==
                                rhs.num_coordinates() * triple_2.row() + triple_1.row()) {
                            if (((sym.inversion == EVEN) ? -parityL : parityL) !=
                                ((sym.reflection == EVEN) ? parityL * parityJ * parityM
                                                          : -parityL * parityJ * parityM)) {
                                continue; // the parity under inversion and reflection is different
                            }
                            skip_reflection =
                                true; // the parity under inversion and reflection is the same
                        }

                        // In case of permutation and reflection symmetry: check whether the
                        // permutation symmetric state is already reflection symmetric
                        if (sym.permutation != NA && col_1 != col_2 && sym.reflection != NA &&
                            mapping[row] ==
                                rhs.num_coordinates() * triple_2.row() + triple_1.row()) {
                            if (((sym.permutation == EVEN) ? -1 : 1) !=
                                ((sym.reflection == EVEN) ? parityL * parityJ * parityM
                                                          : -parityL * parityJ * parityM)) {
                                continue; // the parity under permutation and reflection is
                                          // different
                            }
                            skip_reflection =
                                true; // the parity under permutation and reflection is the same
                        }

                        // In case of inversion and permutation symmetry: the inversion symmetric
                        // state is already permutation symmetric
                        bool skip_permutation = false;
                        if (sym.inversion != NA && sym.permutation != NA && col_1 != col_2) {
                            if (((sym.inversion == EVEN) ? -parityL : parityL) !=
                                ((sym.permutation == EVEN) ? -1 : 1)) {
                                continue; // the parity under inversion and permutation is different
                            }
                            skip_permutation =
                                true; // the parity under inversion and permutation is the same
                        }

                        // In case of rotation symmetry: skip coordinates with wrong total magnetic
                        // momentum
                        if (sym.rotation != NA && sym.rotation != M &&
                            (sym.reflection == NA || sym.rotation != -M)) {
                            continue;
                        }

                        // In case of reflection symmetry: skip half of the coordinates
                        if (sym.reflection != NA && state.m[0] < 0 && !skip_reflection) {
                            continue;
                        }

                        // Calculate coefficient that belongs to the current coordinate
                        scalar_t val_basis = triple_1.value() * triple_2.value(); // coefficient
                        if (sym.reflection != NA && !skip_reflection) {
                            val_basis /= std::sqrt(2);
                        }
                        if (sym.inversion != NA && col_1 != col_2) {
                            val_basis /= std::sqrt(2);
                        }
                        if (sym.permutation != NA && col_1 != col_2 && !skip_permutation) {
                            val_basis /= std::sqrt(2);
                        }

                        // Save the coefficient taking into account the symmetrization
                        mat.addBasis(row, col, val_basis);

                        if (sym.reflection != NA && !skip_reflection) {
                            size_t r = mapping[row];
                            scalar_t v = val_basis;
                            v *= (sym.reflection == EVEN) ? parityL * parityJ * parityM
                                                          : -parityL * parityJ * parityM;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.inversion != NA && col_1 != col_2) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            scalar_t v = val_basis;
                            v *= (sym.inversion == EVEN) ? -parityL : parityL;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.inversion != NA && col_1 != col_2 && sym.reflection != NA &&
                            !skip_reflection) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            r = mapping[r];
                            scalar_t v = val_basis;
                            v *= (sym.reflection == EVEN) ? parityL * parityJ * parityM
                                                          : -parityL * parityJ * parityM;
                            v *= (sym.inversion == EVEN) ? -parityL : parityL;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.permutation != NA && col_1 != col_2 && !skip_permutation) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            scalar_t v = val_basis;
                            v *= (sym.permutation == EVEN) ? -1 : 1;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.permutation != NA && col_1 != col_2 && !skip_permutation &&
                            sym.reflection != NA && !skip_reflection) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            r = mapping[r];
                            scalar_t v = val_basis;
                            v *= (sym.reflection == EVEN) ? parityL * parityJ * parityM
                                                          : -parityL * parityJ * parityM;
                            v *= (sym.permutation == EVEN) ? -1 : 1;
                            mat.addBasis(r, col, v);
                        }

                        existing = true;
                    }
                }

                // Save the combined diagonal element if the corresponding combined basis vector is
                // valid
                if (existing) {
                    mat.addEntries(col, col, val_entries);
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

void energycutoff(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const double &deltaE,
                  std::vector<bool> &necessary) {
    eigen_vector_t diag1 = lhs.entries().diagonal();
    eigen_vector_t diag2 = rhs.entries().diagonal();

    for (int col_1 = 0; col_1 < lhs.basis().outerSize();
         ++col_1) { // outerSize() == num_cols = num_basisvectors()
        for (int col_2 = 0; col_2 < rhs.basis().outerSize(); ++col_2) {
            scalar_t val_entries = diag1[col_1] + diag2[col_2]; // diag(V) x I + I x diag(V)
            if (std::abs(val_entries) < deltaE + 1e-11 ||
                deltaE < 0) { // TODO make +1e-11 unnecessary
                for (eigen_iterator_t triple_1(lhs.basis(), col_1); triple_1; ++triple_1) {
                    for (eigen_iterator_t triple_2(rhs.basis(), col_2); triple_2; ++triple_2) {
                        size_t row =
                            rhs.num_coordinates() * triple_1.row() + triple_2.row(); // coordinate
                        necessary[row] = true;
                    }
                }
            }
        }
    }
}
