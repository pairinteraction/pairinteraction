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

#include "Hamiltonianmatrix.hpp"

#include "EigenCompat.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <fmt/format.h>

#include <stdexcept>

template <typename Scalar>
Hamiltonianmatrix<Scalar>::Hamiltonianmatrix() = default;

template <typename Scalar>
Hamiltonianmatrix<Scalar>::Hamiltonianmatrix(const Eigen::SparseMatrix<Scalar> &entries,
                                             const Eigen::SparseMatrix<Scalar> &basis)
    : entries_(entries), basis_(basis) {}

template <typename Scalar>
Hamiltonianmatrix<Scalar>::Hamiltonianmatrix(size_t szBasis, size_t szEntries) {
    triplets_basis.reserve(szBasis);
    triplets_entries.reserve(szEntries);
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar> &Hamiltonianmatrix<Scalar>::entries() {
    bytes.clear();
    return entries_;
}

template <typename Scalar>
const Eigen::SparseMatrix<Scalar> &Hamiltonianmatrix<Scalar>::entries() const {
    return entries_;
}

template <typename Scalar>
Eigen::SparseMatrix<Scalar> &Hamiltonianmatrix<Scalar>::basis() {
    bytes.clear();
    return basis_;
}

template <typename Scalar>
const Eigen::SparseMatrix<Scalar> &Hamiltonianmatrix<Scalar>::basis() const {
    return basis_;
}

template <typename Scalar>
size_t Hamiltonianmatrix<Scalar>::num_basisvectors() const {
    return basis_.cols();
}

template <typename Scalar>
size_t Hamiltonianmatrix<Scalar>::num_coordinates() const {
    return basis_.rows();
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::addBasis(idx_t row, idx_t col, Scalar val) {
    triplets_basis.emplace_back(row, col, val);
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::addEntries(idx_t row, idx_t col, Scalar val) {
    triplets_entries.emplace_back(row, col, val);
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::compress(size_t nBasis, size_t nCoordinates) {
    basis_.resize(nCoordinates, nBasis);
    entries_.resize(nBasis, nBasis);
    basis_.setFromTriplets(triplets_basis.begin(), triplets_basis.end());
    entries_.setFromTriplets(triplets_entries.begin(), triplets_entries.end());
    triplets_basis.clear();
    triplets_entries.clear();
}

template <typename Scalar>
std::vector<Hamiltonianmatrix<Scalar>> Hamiltonianmatrix<Scalar>::findSubs() const { // TODO
    std::vector<Hamiltonianmatrix<Scalar>> submatrices;
    submatrices.push_back(*this);
    return submatrices;
}

template <typename Scalar>
Hamiltonianmatrix<Scalar> Hamiltonianmatrix<Scalar>::abs() const {
    return Hamiltonianmatrix<Scalar>(entries_.cwiseAbs().template cast<Scalar>(), basis_);
}

template <typename Scalar>
Hamiltonianmatrix<Scalar>
Hamiltonianmatrix<Scalar>::changeBasis(const Eigen::SparseMatrix<Scalar> &basis) const {
    auto transformator = basis_.adjoint() * basis;
    auto entries = transformator.adjoint() * entries_ * transformator;
    return Hamiltonianmatrix<Scalar>(entries, basis);
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::applyCutoff(double cutoff) {
    bytes.clear();

    // build transformator
    Eigen::VectorX<Scalar> diag = entries_.diagonal();

    std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
    triplets_transformator.reserve(num_basisvectors());

    size_t idxBasis = 0;
    for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
        if (std::abs(diag[idx]) < cutoff) {
            triplets_transformator.emplace_back(idx, idxBasis++, 1);
        }
    }

    Eigen::SparseMatrix<Scalar> transformator(this->num_basisvectors(), idxBasis);
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = basis_ * transformator;
    entries_ = transformator.adjoint() * entries_ * transformator;
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::findUnnecessaryStates(
    std::vector<bool> &isNecessaryCoordinate) const {
    std::vector<double> isNecessaryCoordinate_real(num_coordinates(), 0);
    for (int k = 0; k < basis_.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basis_, k); triple;
             ++triple) {
            isNecessaryCoordinate_real[triple.row()] += std::pow(std::abs(triple.value()), 2);
        }
    }

    for (size_t idx = 0; idx < this->num_coordinates(); ++idx) {
        if (isNecessaryCoordinate_real[idx] > 0.05) { // TODO
            isNecessaryCoordinate[idx] = true;
        }
    }
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::removeUnnecessaryBasisvectors(
    const std::vector<bool> &isNecessaryCoordinate) {
    bytes.clear();

    // build transformator
    std::vector<double> isNecessaryBasisvector(num_basisvectors(), 0);
    for (int k_1 = 0; k_1 < basis_.outerSize(); ++k_1) {
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basis_, k_1); triple;
             ++triple) {
            ptrdiff_t col = triple.col(); // basis vector
            ptrdiff_t row = triple.row(); // coordinate
            if (isNecessaryCoordinate[row]) {
                isNecessaryBasisvector[col] += std::pow(std::abs(triple.value()), 2);
            }
        }
    }

    std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
    triplets_transformator.reserve(num_basisvectors());

    size_t idxBasis = 0;
    for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
        if (isNecessaryBasisvector[idx] > 0.05) { // TODO
            triplets_transformator.emplace_back(idx, idxBasis++, 1);
        }
    }

    Eigen::SparseMatrix<Scalar> transformator(this->num_basisvectors(), idxBasis);
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = basis_ * transformator;
    entries_ = transformator.adjoint() * entries_ * transformator;
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::removeUnnecessaryBasisvectors() {
    bytes.clear();

    // build transformator
    std::vector<double> isNecessaryBasisvector(num_basisvectors(), 0);
    for (int k_1 = 0; k_1 < basis_.outerSize(); ++k_1) {
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple(basis_, k_1); triple;
             ++triple) {
            ptrdiff_t col = triple.col(); // basis vector
            isNecessaryBasisvector[col] += std::pow(std::abs(triple.value()), 2);
        }
    }

    std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
    triplets_transformator.reserve(num_basisvectors());

    size_t idxBasis = 0;
    for (size_t idx = 0; idx < this->num_basisvectors(); ++idx) {
        if (isNecessaryBasisvector[idx] > 0.05) {
            triplets_transformator.emplace_back(idx, idxBasis++, 1);
        }
    }

    Eigen::SparseMatrix<Scalar> transformator(this->num_basisvectors(), idxBasis);
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = basis_ * transformator;
    entries_ = transformator.adjoint() * entries_ * transformator;
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::removeUnnecessaryStates(
    const std::vector<bool> &isNecessaryCoordinate) {
    bytes.clear();

    // build transformator
    std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
    triplets_transformator.reserve(num_coordinates());

    size_t idxCoordinate = 0;
    for (size_t idx = 0; idx < this->num_coordinates(); ++idx) {
        if (isNecessaryCoordinate[idx]) {
            triplets_transformator.emplace_back(idxCoordinate++, idx, 1);
        }
    }

    Eigen::SparseMatrix<Scalar> transformator(idxCoordinate, this->num_coordinates());
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    // apply transformator
    basis_ = transformator * basis_;
}

template <typename Scalar>
Hamiltonianmatrix<Scalar>
Hamiltonianmatrix<Scalar>::getBlock(const std::vector<ptrdiff_t> &indices) {
    std::vector<Eigen::Triplet<Scalar>> triplets_transformator;
    triplets_transformator.reserve(indices.size());
    for (size_t idx = 0; idx < indices.size(); ++idx) {
        triplets_transformator.emplace_back(indices[idx], idx, 1);
    }
    Eigen::SparseMatrix<Scalar> transformator(this->num_basisvectors(), indices.size());
    transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

    Eigen::SparseMatrix<Scalar> block_entries = transformator.adjoint() * entries_ * transformator;
    Eigen::SparseMatrix<Scalar> block_basis = basis_ * transformator;

    return Hamiltonianmatrix<Scalar>(block_entries, block_basis);
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::diagonalize() {
    if (this->num_basisvectors() > 1) { // NOLINT
        // diagonalization
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixX<Scalar>> eigensolver(
            Eigen::MatrixX<Scalar>(this->entries()));

        // eigenvalues and eigenvectors
        Eigen::VectorX<double> evals = eigensolver.eigenvalues();
        Eigen::SparseMatrix<Scalar> evecs = eigensolver.eigenvectors().sparseView(1e-4, 0.5);

        this->entries().setZero();
        this->entries().reserve(evals.size());
        for (int idx = 0; idx < evals.size(); ++idx) {
            this->entries().insert(idx, idx) = evals.coeffRef(idx);
        }
        this->entries().makeCompressed();

        this->basis() = (this->basis() * evecs).pruned(1e-4, 0.5);
    }
}

template <typename Scalar>
Hamiltonianmatrix<Scalar> operator+(Hamiltonianmatrix<Scalar> lhs,
                                    const Hamiltonianmatrix<Scalar> &rhs) {
    lhs.bytes.clear();
    lhs.entries_ += rhs.entries_;
    return lhs;
}

template <typename Scalar>
Hamiltonianmatrix<Scalar> operator-(Hamiltonianmatrix<Scalar> lhs,
                                    const Hamiltonianmatrix<Scalar> &rhs) {
    lhs.bytes.clear();
    lhs.entries_ -= rhs.entries_;
    return lhs;
}

template <typename Scalar, typename T>
Hamiltonianmatrix<Scalar> operator*(const T &lhs, Hamiltonianmatrix<Scalar> rhs) {
    rhs.bytes.clear();
    rhs.entries_ *= lhs;
    return rhs;
}

template <typename Scalar, typename T>
Hamiltonianmatrix<Scalar> operator*(Hamiltonianmatrix<Scalar> lhs, const T &rhs) {
    lhs.bytes.clear();
    lhs.entries_ *= rhs;
    return lhs;
}

template <typename Scalar>
Hamiltonianmatrix<Scalar> &
Hamiltonianmatrix<Scalar>::operator+=(const Hamiltonianmatrix<Scalar> &rhs) {
    bytes.clear();
    entries_ += rhs.entries_;
    return *this;
}

template <typename Scalar>
Hamiltonianmatrix<Scalar> &
Hamiltonianmatrix<Scalar>::operator-=(const Hamiltonianmatrix<Scalar> &rhs) {
    bytes.clear();
    entries_ -= rhs.entries_;
    return *this;
}

template <typename Scalar>
bytes_t &Hamiltonianmatrix<Scalar>::serialize() {
    doSerialization();
    return bytes;
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::doSerialization() {
    if (bytes.empty()) {
        entries_.makeCompressed();
        basis_.makeCompressed();

        // convert matrix "entries" to vectors of primitive data types
        byte_t entries_flags = 0;
        if (entries_.IsRowMajor != 0) {
            entries_flags |= csr_not_csc;
        }
        if (utils::is_complex<Scalar>::value) {
            entries_flags |= complex_not_real;
        }
        storage_idx_t entries_rows = entries_.rows();
        storage_idx_t entries_cols = entries_.cols();
        std::vector<Scalar> entries_data(entries_.valuePtr(),
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
        if (utils::is_complex<Scalar>::value) {
            basis_flags |= complex_not_real;
        }
        storage_idx_t basis_rows = basis_.rows();
        storage_idx_t basis_cols = basis_.cols();
        std::vector<Scalar> basis_data(basis_.valuePtr(), basis_.valuePtr() + basis_.nonZeros());
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

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::deserialize(bytes_t &bytesin) {
    bytes = bytesin;
    doDeserialization();
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::doDeserialization() {
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

    if ((((entries_flags & complex_not_real) > 0) != utils::is_complex<Scalar>::value) ||
        (((basis_flags & complex_not_real) > 0) != utils::is_complex<Scalar>::value)) {
        std::string msg("The data type used in the program does not fit the data type used in the "
                        "serialized objects.");
        std::cout << fmt::format(">>ERR{:s}", msg.c_str()) << std::endl;
        throw std::runtime_error(msg);
    }

    // build matrix "entries_"
    std::vector<Scalar> entries_data;
    mergeComplex(entries_data_real, entries_data_imag, entries_data);
    entries_ = Eigen::SparseMatrix<Scalar>(entries_rows, entries_cols);
    entries_.makeCompressed();
    entries_.resizeNonZeros(entries_data.size());
    std::copy(entries_data.begin(), entries_data.end(), entries_.valuePtr());
    std::copy(entries_indices.begin(), entries_indices.end(), entries_.innerIndexPtr());
    std::copy(entries_indptr.begin(), entries_indptr.end(), entries_.outerIndexPtr());
    entries_.finalize();

    // build matrix "basis_"
    std::vector<Scalar> basis_data;
    mergeComplex(basis_data_real, basis_data_imag, basis_data);
    basis_ = Eigen::SparseMatrix<Scalar>(basis_rows, basis_cols);
    basis_.makeCompressed();
    basis_.resizeNonZeros(basis_data.size());
    std::copy(basis_data.begin(), basis_data.end(), basis_.valuePtr());
    std::copy(basis_indices.begin(), basis_indices.end(), basis_.innerIndexPtr());
    std::copy(basis_indptr.begin(), basis_indptr.end(), basis_.outerIndexPtr());
    basis_.finalize();
}

template <typename Scalar>
uint64_t Hamiltonianmatrix<Scalar>::hashEntries() {
    // TODO bring this functionality to the matrix class and use it for serialization, too
    doSerialization();
    return utils::FNV64(&bytes[0], bytes.size());
}

template <typename Scalar>
uint64_t Hamiltonianmatrix<Scalar>::hashBasis() {
    // TODO bring this functionality to the matrix class and use it for serialization, too
    doSerialization();
    return utils::FNV64(&bytes[0], bytes.size());
}

template <typename Scalar>
void Hamiltonianmatrix<Scalar>::save(const std::string &fname) {
    doSerialization();

    // open file
    FILE *pFile;
    pFile = fopen(fname.c_str(), "wb");

    // write
    fwrite(&bytes[0], 1, sizeof(byte_t) * bytes.size(), pFile);

    // close file
    fclose(pFile);
}

template <typename Scalar>
bool Hamiltonianmatrix<Scalar>::load(const std::string &fname) {
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

template <typename Scalar>
Hamiltonianmatrix<Scalar> combine(const Hamiltonianmatrix<Scalar> &lhs,
                                  const Hamiltonianmatrix<Scalar> &rhs, const double &deltaE,
                                  const std::shared_ptr<BasisnamesTwo> &basis_two,
                                  const Symmetry &sym) {
    // TODO program a faster method for samebasis == true

    size_t num_basisvectors = lhs.num_basisvectors() * rhs.num_basisvectors();
    size_t num_coordinates = lhs.num_coordinates() * rhs.num_coordinates();

    ////////////////////////////////////////////////////////
    ////// Mapping used in case of reflection symmetry /////
    ////////////////////////////////////////////////////////

    std::vector<size_t> mapping(num_coordinates, -1);
    if (sym.reflection != NA) { // NOLINT
        std::unordered_map<StateTwoOld, size_t> buffer;
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
    Eigen::VectorX<Scalar> diag1 = lhs.entries().diagonal();
    Eigen::VectorX<Scalar> diag2 = rhs.entries().diagonal();

    // Number of elements for which space sould be reserved
    size_t size_basis = num_basisvectors; // TODO estimate better
    size_t size_entries = num_basisvectors;

    Hamiltonianmatrix<Scalar> mat(size_basis, size_entries);

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
            Scalar val_entries = diag1[col_1] + diag2[col_2]; // diag(V) x I + I x diag(V)

            // Check whether the new diagonal element is within the energy cutoff
            if (std::abs(val_entries) < deltaE + 1e-11 ||
                deltaE < 0) { // TODO avoid the "+1e-11" hack

                // Variable that stores whether the combindes basis vector, that belongs to the
                // combined diagonal element, is valid
                bool existing = false;

                // --- Combine basis vectors for mat.basis() ---
                for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple_1(lhs.basis(),
                                                                                  col_1);
                     triple_1; ++triple_1) {
                    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple_2(rhs.basis(),
                                                                                      col_2);
                         triple_2; ++triple_2) {
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
                        Scalar val_basis = triple_1.value() * triple_2.value(); // coefficient
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
                            Scalar v = val_basis;
                            v *= (sym.reflection == EVEN) ? parityL * parityJ * parityM
                                                          : -parityL * parityJ * parityM;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.inversion != NA && col_1 != col_2) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            Scalar v = val_basis;
                            v *= (sym.inversion == EVEN) ? -parityL : parityL;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.inversion != NA && col_1 != col_2 && sym.reflection != NA &&
                            !skip_reflection) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            r = mapping[r];
                            Scalar v = val_basis;
                            v *= (sym.reflection == EVEN) ? parityL * parityJ * parityM
                                                          : -parityL * parityJ * parityM;
                            v *= (sym.inversion == EVEN) ? -parityL : parityL;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.permutation != NA && col_1 != col_2 && !skip_permutation) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            Scalar v = val_basis;
                            v *= (sym.permutation == EVEN) ? -1 : 1;
                            mat.addBasis(r, col, v);
                        }

                        if (sym.permutation != NA && col_1 != col_2 && !skip_permutation &&
                            sym.reflection != NA && !skip_reflection) {
                            size_t r = rhs.num_coordinates() * triple_2.row() + triple_1.row();
                            r = mapping[r];
                            Scalar v = val_basis;
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

template <typename Scalar>
void energycutoff(const Hamiltonianmatrix<Scalar> &lhs, const Hamiltonianmatrix<Scalar> &rhs,
                  const double &deltaE, std::vector<bool> &necessary) {
    Eigen::VectorX<Scalar> diag1 = lhs.entries().diagonal();
    Eigen::VectorX<Scalar> diag2 = rhs.entries().diagonal();

    for (int col_1 = 0; col_1 < lhs.basis().outerSize();
         ++col_1) { // outerSize() == num_cols = num_basisvectors()
        for (int col_2 = 0; col_2 < rhs.basis().outerSize(); ++col_2) {
            Scalar val_entries = diag1[col_1] + diag2[col_2]; // diag(V) x I + I x diag(V)
            if (std::abs(val_entries) < deltaE + 1e-11 ||
                deltaE < 0) { // TODO make +1e-11 unnecessary
                for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple_1(lhs.basis(),
                                                                                  col_1);
                     triple_1; ++triple_1) {
                    for (typename Eigen::SparseMatrix<Scalar>::InnerIterator triple_2(rhs.basis(),
                                                                                      col_2);
                         triple_2; ++triple_2) {
                        size_t row =
                            rhs.num_coordinates() * triple_1.row() + triple_2.row(); // coordinate
                        necessary[row] = true;
                    }
                }
            }
        }
    }
}

template class Hamiltonianmatrix<std::complex<double>>;
template Hamiltonianmatrix<std::complex<double>>
operator+(Hamiltonianmatrix<std::complex<double>> lhs,
          const Hamiltonianmatrix<std::complex<double>> &rhs);
template Hamiltonianmatrix<std::complex<double>>
operator-(Hamiltonianmatrix<std::complex<double>> lhs,
          const Hamiltonianmatrix<std::complex<double>> &rhs);
template Hamiltonianmatrix<std::complex<double>>
operator*(const std::complex<double> &lhs, Hamiltonianmatrix<std::complex<double>> rhs);
template Hamiltonianmatrix<std::complex<double>>
operator*(Hamiltonianmatrix<std::complex<double>> lhs, const std::complex<double> &rhs);
template Hamiltonianmatrix<std::complex<double>>
operator*(const double &lhs, Hamiltonianmatrix<std::complex<double>> rhs);
template Hamiltonianmatrix<std::complex<double>>
operator*(Hamiltonianmatrix<std::complex<double>> lhs, const double &rhs);
template Hamiltonianmatrix<std::complex<double>>
combine(const Hamiltonianmatrix<std::complex<double>> &lhs,
        const Hamiltonianmatrix<std::complex<double>> &rhs, const double &deltaE,
        const std::shared_ptr<BasisnamesTwo> &basis_two, const Symmetry &sym);
template void energycutoff(const Hamiltonianmatrix<std::complex<double>> &lhs,
                           const Hamiltonianmatrix<std::complex<double>> &rhs, const double &deltaE,
                           std::vector<bool> &necessary);
template class Hamiltonianmatrix<double>;
template Hamiltonianmatrix<double> operator+(Hamiltonianmatrix<double> lhs,
                                             const Hamiltonianmatrix<double> &rhs);
template Hamiltonianmatrix<double> operator-(Hamiltonianmatrix<double> lhs,
                                             const Hamiltonianmatrix<double> &rhs);
template Hamiltonianmatrix<double> operator*(const double &lhs, Hamiltonianmatrix<double> rhs);
template Hamiltonianmatrix<double> operator*(Hamiltonianmatrix<double> lhs, const double &rhs);
template Hamiltonianmatrix<double>
combine(const Hamiltonianmatrix<double> &lhs, const Hamiltonianmatrix<double> &rhs,
        const double &deltaE, const std::shared_ptr<BasisnamesTwo> &basis_two, const Symmetry &sym);
template void energycutoff(const Hamiltonianmatrix<double> &lhs,
                           const Hamiltonianmatrix<double> &rhs, const double &deltaE,
                           std::vector<bool> &necessary);
