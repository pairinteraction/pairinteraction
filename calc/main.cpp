#include "dtypes.h"
#include "MpiEnvironment.h"
#include "MpiLoadbalancingComplex.h"
#include "MpiLoadbalancingSimple.h"
#include "Vectorizable.h"
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
#include <getopt.h>
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

#include <numeric>


/*
///////////////////// Sources /////////////////////

* http://cpplove.blogspot.co.at/2013/05/my-take-on-serialization-part-ii.html
* http://www.codeproject.com/Articles/567242/AplusC-b-bplusObjectplusFactory
* https://stackoverflow.com/questions/5120768/how-to-implement-the-factory-method-pattern-in-c-correctly


///////////////////// TODOs /////////////////////

* use Google C++ Style Guide (http://google.github.io/styleguide/cppguide.html)
* mpi namespace
* #define MPI_FLAG    #ifdef MPI_FLAG ... #endif
* use std::vector< bytes_t> instead of std::vector<triple>
* do not use the mpi c++ binding (it is depreciated)

* second loop only over relevant basis states
* find submatrices
* parallelize construction of Hamiltonian
* nur eine Hälfte der symmetrischen Matrizen aufbauen


///////////////////// Structure /////////////////////

// Get start state from json file
State startstate(json);

// Decompose start state into non-interacting states (using asym/sym, Wigner d-matrices)
superposition = startstate.decompose();

// Figure out whether the states of the superposition are already fully cached
superposition.getCached(vecR, vecB, vecE);

// Create basis
Basises basises(superposition); // use loadbalancing without (*) slave groups, do nothing for fully cached states

// Create Hamiltonian
Hamiltonians hamiltonians(basises, superposition); // use loadbalancing without (*) slave groups, extract relevant submatrices, do nothing for fully cached states

// Create Eigensystem
Eigensystems eigensystems(hamiltonians, vecR, vecB, vecE); // use loadbalancing WITH slave groups, read and write to cache if necessary
eigensystems.save(json);

// Analyze
potentials = eigensystems.calculatePotentials();
potentials.save(json);

overlap = eigensystems.calculateOverlap(superposition);
overlap.save(json);

(*) maybe WITH slave groups, too


///////////////////// Commands /////////////////////

char filename[20+1+3+1]; sprintf(filename, "%020" PRIu64 ".mat", FNV64(&bytes[idxStart], bytes.size()));
*/




// ----------------------------------------

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
        auto transformator = basis_.transpose()*basis; // TODO conjugate transpose / inverse (auch an anderen Stellen) !!!!!!!!!!!!!!!!!
        auto entries = transformator.transpose()*entries_*transformator;
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
        entries_= transformator.transpose()*entries_*transformator;
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
        entries_= transformator.transpose()*entries_*transformator;
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
            if (isNecessaryBasisvector[idx] > 0.05) { // TODO
                triplets_transformator.push_back(eigen_triplet_t(idx,idxBasis++,1));
            }
        }

        eigen_sparse_t transformator(this->num_basisvectors(),idxBasis);
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // apply transformator
        basis_ = basis_*transformator;
        entries_= transformator.transpose()*entries_*transformator;
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

        eigen_sparse_t block_entries = transformator.transpose()*entries_*transformator;
        eigen_sparse_t block_basis = basis_*transformator;

        return Hamiltonianmatrix(block_entries,  block_basis);
    }

    friend Hamiltonianmatrix combineSym(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const real_t &deltaE, const std::vector<bool> &necessary1, const std::vector<bool> &necessary2, const std::vector<bool> &necessary) {
        return lhs.duplicate(SYM, rhs, deltaE, necessary1, necessary2, necessary);
    }
    friend Hamiltonianmatrix combineAsym(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const real_t &deltaE, const std::vector<bool> &necessary1, const std::vector<bool> &necessary2, const std::vector<bool> &necessary) {
        return lhs.duplicate(ASYM, rhs, deltaE, necessary1, necessary2, necessary);
    }
    friend Hamiltonianmatrix combineAll(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs, const real_t &deltaE, const std::vector<bool> &necessary1, const std::vector<bool> &necessary2, const std::vector<bool> &necessary) {
        return lhs.duplicate(ALL, rhs, deltaE, necessary1, necessary2, necessary);
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

            // convert filename to vector of primitive data type
            std::vector<char> name(filename_.begin(), filename_.end());

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
            s << name;
            s << isExisting_;
            idxStart = s.position();
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
        std::vector<char> name;
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
        s >> name;
        s >> isExisting_;
        idxStart = s.position();
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
            std::cout << ">>ERR" << "The data type used in the program does not fit the data type used in the serialized objects." << std::endl;
            abort();
        }

        // build filename
        filename_ = std::string(&name[0], name.size());

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
    void addFilename(std::string name) {
        bytes.clear();
        filename_ = name;
    }
    void addIsExisting(bool isExisting) {
        bytes.clear();
        isExisting_ = isExisting;
    }

    std::string& filename() {
        return filename_;
    }

    uint64_t hashEntries() {
        // TODO bring this functionality to the matrix class and use it for serialization, too
        doSerialization();
        return utils::FNV64(&bytes[idxStart], bytes.size());
    }

    uint64_t hashBasis() {
        // TODO bring this functionality to the matrix class and use it for serialization, too
        doSerialization();
        return utils::FNV64(&bytes[idxStart], bytes.size());
    }

    bool exist() {
        return isExisting_;
    }

    void save() {
        doSerialization();

        // open file
        FILE *pFile;
        pFile = fopen(filename_.c_str() , "wb" );

        // write
        fwrite(&bytes[idxStart], 1 , sizeof(byte_t)*(bytes.size()-idxStart), pFile );

        // close file
        fclose(pFile);
    }
    bool load() {
        // open file
        if (FILE *pFile = fopen(filename_.c_str() , "rb" )) {

            // obtain file size:
            fseek (pFile , 0 , SEEK_END);
            size_t size_file = ftell(pFile);
            rewind(pFile);

            // prepare bytes
            std::vector<char> name(filename_.begin(), filename_.end());
            Serializer s;
            s << name;
            s << isExisting_;
            idxStart = s.position();
            s.save(bytes);
            bytes.resize(size_file/sizeof(byte_t)+idxStart);

            // read
            size_t size_result = fread(&bytes[idxStart], 1 , sizeof(byte_t)*(bytes.size()-idxStart) , pFile );
            if (size_result != size_file) throw std::runtime_error("Matrix could not be read from file.");

            // close file
            fclose(pFile);

            doDeserialization();
            return true;
        } else {
            return false;
        }
    }

protected:
    enum mode_t {ALL, SYM, ASYM};
    void addSymmetrized(mode_t mode, std::vector<idx_t> mapping, idx_t row_1, idx_t row_2, idx_t col_1, idx_t col_2, scalar_t val, std::vector<eigen_triplet_t> &triplets_entries) const {
        idx_t row = mapping[this->num_basisvectors()*row_1 + row_2];
        idx_t col = mapping[this->num_basisvectors()*col_1 + col_2];
        if((mode == ALL) || (mode == SYM && row_1 <= row_2 && col_1 <= col_2) || (mode == ASYM && row_1 < row_2 && col_1 < col_2)) {
            real_t factor = 1;
            if (mode == SYM && row_1 == row_2) factor *= 1./sqrt(2.);
            if (mode == SYM && col_1 == col_2) factor *= 1./sqrt(2.);
            triplets_entries.push_back(eigen_triplet_t(row, col, factor*val));
        }
    }
    Hamiltonianmatrix duplicate(mode_t mode, const Hamiltonianmatrix &rhs, const real_t &deltaE, const std::vector<bool> &necessary1, const std::vector<bool> &necessary2, const std::vector<bool> &necessary) const {
        real_t tol = 1e-32;

        size_t num_basisvectors = this->num_basisvectors()*rhs.num_basisvectors();
        size_t num_coordinates = this->num_coordinates()*rhs.num_coordinates();

        // --- mapping ---
        std::vector<ptrdiff_t> mapping(num_basisvectors, -1);

        eigen_vector_t diag1 = entries_.diagonal();
        eigen_vector_t diag2 = rhs.entries().diagonal();

        size_t i = 0;

        for (size_t idx_1 = 0; idx_1 < this->num_basisvectors(); ++idx_1) {
            for (size_t idx_2 = 0; idx_2 < rhs.num_basisvectors(); ++idx_2) {
                if (mode == SYM && idx_2 < idx_1) {
                    continue;
                }

                if (mode == ASYM && idx_2 <= idx_1) {
                    continue;
                }

                size_t idx = rhs.num_basisvectors()*idx_1 + idx_2;
                scalar_t val = diag1[idx_1] + diag2[idx_2]; // diag(V) x I + I x diag(V)

                if (std::abs(val) < deltaE+1e-11 || deltaE < 0) {
                    mapping[idx] = i++;
                }
            }
        }

        size_t num_basisvectors_new = i;

        // --- initialize matrix ---
        size_t size_basis = basis_.nonZeros() / (4.*num_basisvectors) * rhs.basis().nonZeros() * num_basisvectors_new ; //basis_.nonZeros()/this->num_basisvectors() * rhs.basis().nonZeros()/rhs.num_basisvectors() * num_basisvectors;
        size_t size_entries = num_basisvectors;

        Hamiltonianmatrix mat(size_basis, size_entries);

        num_basisvectors = num_basisvectors_new;

        // --- duplicate basis_ ---

        for (eigen_idx_t k_1=0; k_1<basis_.outerSize(); ++k_1) {

            for (eigen_iterator_t triple_1(basis_,k_1); triple_1; ++triple_1) {
                if (!necessary1[triple_1.row()]) continue;

                for (eigen_idx_t k_2=0; k_2<rhs.basis().outerSize(); ++k_2) {

                    for (eigen_iterator_t triple_2(rhs.basis(),k_2); triple_2; ++triple_2) {
                        if (!necessary2[triple_2.row()]) continue;

                        ptrdiff_t col = mapping[rhs.num_basisvectors()*triple_1.col() + triple_2.col()]; // basis vector

                        if (col >= 0) {
                            if (mode == ALL || triple_1.col() == triple_2.col()) {
                                size_t row = rhs.num_coordinates()*triple_1.row() + triple_2.row(); // coordinate
                                scalar_t val = triple_1.value() * triple_2.value();
                                if (necessary[row]) mat.addBasis(row,col,val);

                            } else {
                                size_t row = rhs.num_coordinates()*triple_1.row() + triple_2.row(); // coordinate
                                scalar_t val = triple_1.value() * triple_2.value();
                                val /= std::sqrt(2);
                                if (necessary[row]) mat.addBasis(row,col,val);

                                row = rhs.num_coordinates()*triple_2.row() + triple_1.row(); // coordinate
                                val *= (mode == ASYM) ? -1 : 1;
                                if (necessary[row]) mat.addBasis(row,col,val);
                            }

                        }
                    }
                }
            }
        }

        // --- duplicate entries_ (does proberly work only if entries_ is diagonal) --- // TODO

        // V x I
        for (eigen_idx_t k=0; k<entries_.outerSize(); ++k) {
            for (eigen_iterator_t it(entries_,k); it; ++it) {
                for (size_t unitmatrix_idx = 0; unitmatrix_idx < rhs.num_basisvectors(); ++unitmatrix_idx) {
                    scalar_t val = it.value();

                    if (std::abs(val) > tol) {
                        size_t row_1 = it.row();
                        size_t row_2 = unitmatrix_idx;
                        size_t col_1 = it.col();
                        size_t col_2 = unitmatrix_idx;

                        ptrdiff_t row = mapping[rhs.num_basisvectors()*row_1 + row_2];
                        ptrdiff_t col = mapping[rhs.num_basisvectors()*col_1 + col_2];

                        if (row >= 0 && col >= 0) {
                            mat.addEntries(row,col,val);
                        }
                    }
                }
            }
        }

        // I x V
        for (eigen_idx_t k=0; k<rhs.entries().outerSize(); ++k) {
            for (eigen_iterator_t it(rhs.entries(),k); it; ++it) {
                for (size_t unitmatrix_idx = 0; unitmatrix_idx < this->num_basisvectors(); ++unitmatrix_idx) {
                    scalar_t val = it.value();

                    if (std::abs(val) > tol) {
                        size_t row_1 = unitmatrix_idx;
                        size_t row_2 = it.row();
                        size_t col_1 = unitmatrix_idx;
                        size_t col_2 = it.col();

                        ptrdiff_t row = mapping[rhs.num_basisvectors()*row_1 + row_2];
                        ptrdiff_t col = mapping[rhs.num_basisvectors()*col_1 + col_2];

                        if (row >= 0 && col >= 0) {
                            mat.addEntries(row,col,val);
                        }
                    }
                }
            }
        }

        mat.compress(num_basisvectors, num_coordinates);

        return mat;
    }
    Hamiltonianmatrix duplicate(mode_t mode, const Hamiltonianmatrix &rhs) const {
        // --- mapping ---
        idx_t i = 0;
        std::vector<idx_t> mapping(this->num_basisvectors()*rhs.num_basisvectors());
        for (idx_t idx_1 = 0; idx_1 < this->num_basisvectors(); ++idx_1) {
            for (idx_t idx_2 = 0; idx_2 < rhs.num_basisvectors(); ++idx_2) {
                if ((mode != ALL && idx_1 < idx_2) || (mode == ALL) || (mode == SYM && idx_1 == idx_2)) {
                    idx_t idx = rhs.num_coordinates()*idx_1 + idx_2;
                    mapping[idx] = i++;
                }
            }
        }

        Hamiltonianmatrix mat(i,this->num_coordinates()*rhs.num_coordinates());

        // --- duplicate basis_ --- // TODO
        std::vector<eigen_triplet_t> triplets_basis;
        triplets_basis.reserve(basis_.nonZeros()*rhs.basis().nonZeros());

        for (eigen_idx_t k_1=0; k_1<basis_.outerSize(); ++k_1) {
            for (eigen_iterator_t triple_1(basis_,k_1); triple_1; ++triple_1) {
                for (eigen_idx_t k_2=0; k_2<rhs.basis().outerSize(); ++k_2) {
                    for (eigen_iterator_t triple_2(rhs.basis(),k_2); triple_2; ++triple_2) {
                        if ((mode != ALL && triple_1.col() < triple_2.col())) {
                            idx_t idx_row1 = rhs.num_basisvectors()*triple_1.row() + triple_2.row(); // coord1
                            idx_t idx_row2 = this->num_basisvectors()*triple_2.row() + triple_1.row(); // coord2
                            idx_t idx_col = mapping[rhs.num_coordinates()*triple_1.col() + triple_2.col()]; // vec

                            int factor = (mode == ASYM) ? -1 : 1;
                            real_t multiplier = 1/sqrt(2.);
                            triplets_basis.push_back(eigen_triplet_t(idx_row1, idx_col, triple_1.value()*triple_2.value()*multiplier));
                            multiplier = 1/sqrt(2.)*factor;
                            triplets_basis.push_back(eigen_triplet_t(idx_row2, idx_col, triple_1.value()*triple_2.value()*multiplier));

                        } else if ((mode == ALL) || (mode == SYM && triple_1.col() == triple_2.col())) {
                            idx_t idx_row = rhs.num_basisvectors()*triple_1.row() + triple_2.row(); // coord
                            idx_t idx_col = mapping[rhs.num_coordinates()*triple_1.col() + triple_2.col()]; // vec
                            triplets_basis.push_back(eigen_triplet_t(idx_row, idx_col, triple_1.value()*triple_2.value()));
                        }
                    }
                }
            }
        }

        mat.basis().setFromTriplets(triplets_basis.begin(), triplets_basis.end());

        // --- duplicate entries_ --- // TODO
        std::vector<eigen_triplet_t> triplets_entries;
        triplets_entries.reserve(2*entries_.nonZeros()*this->num_basisvectors());

        for (eigen_idx_t k_1=0; k_1<entries_.outerSize(); ++k_1) {
            for (eigen_iterator_t hamiltonian_triple(entries_,k_1); hamiltonian_triple; ++hamiltonian_triple) {
                for (size_t unitmatrix_idx = 0; unitmatrix_idx < this->num_basisvectors(); ++unitmatrix_idx) {
                    idx_t row_1, row_2, col_1, col_2;
                    scalar_t val = hamiltonian_triple.value();

                    // --- ordered terms ---
                    // <1a 2a|V x I|1b 2b> = <2a 1a|I x V|2b 1b>
                    row_1 = hamiltonian_triple.row();
                    row_2 = unitmatrix_idx;
                    col_1 = hamiltonian_triple.col();
                    col_2 = unitmatrix_idx;
                    addSymmetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);

                    // <1a 2a|I x V|1b 2b> = <2a 1a|V x I|2b 1b>
                    row_1 = unitmatrix_idx;
                    row_2 = hamiltonian_triple.row();
                    col_1 = unitmatrix_idx;
                    col_2 = hamiltonian_triple.col();
                    addSymmetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);

                    // --- mixed terms ---
                    if (mode == ALL) continue;

                    // <1a 2a|V x I|2b 1b> = <2a 1a|I x V|1b 2b>
                    row_1 = hamiltonian_triple.row();
                    row_2 = unitmatrix_idx;
                    col_1 = unitmatrix_idx;
                    col_2 = hamiltonian_triple.col();
                    addSymmetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);

                    // <1a 2a|I x V|2b 1b> = <2a 1a|V x I|1b 2b>
                    row_1 = unitmatrix_idx;
                    row_2 = hamiltonian_triple.row();
                    col_1 = hamiltonian_triple.col();
                    col_2 = unitmatrix_idx;
                    addSymmetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);
                }
            }
        }

        mat.entries().setFromTriplets(triplets_entries.begin(), triplets_entries.end());

        return mat;
    }

    eigen_sparse_t entries_;
    eigen_sparse_t basis_;

    bytes_t bytes;
    size_t idxStart;

    std::string filename_;
    bool isExisting_;

    std::vector<eigen_triplet_t> triplets_basis;
    std::vector<eigen_triplet_t> triplets_entries;
};

class Hamiltonian : protected MpiLoadbalancingSimple<Hamiltonianmatrix, Hamiltonianmatrix> {
public:
    Hamiltonian(std::shared_ptr<MpiEnvironment> mpi) : MpiLoadbalancingSimple(mpi, 1000) {}
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
    void saveLines() const {
        if (mpi->rank() == 0) {
            std::shared_ptr<Hamiltonianmatrix> mat_previous = nullptr;

            //std::vector<std::vector<real_t>> lines; // TODO stattdessen indices+1 in sparse matrix abspeichern

            idx_t lines_idx_max = 0;
            std::vector<eigen_triplet_real_t> triplets_lines_eigenenergies;
            triplets_lines_eigenenergies.reserve(matrix_diag.size()*matrix_diag.front()->entries().outerSize()*1.2);
            std::unordered_map<eigen_idx_t,eigen_idx_t> lines_idx;

            eigen_sparse_real_t overlap_root(matrix_diag.front()->basis().rows(),matrix_diag.front()->basis().rows());

            idx_t step = 0;
            for (auto &mat_current: matrix_diag) {
                if (mat_previous != nullptr) {
                    overlap_root = (mat_previous->basis().adjoint() * mat_current->basis()).pruned(1e-4,0.5).cwiseAbs(); // for row major order, the result has to be transposed
                }
                eigen_vector_real_t eigenenergies = mat_current->entries().diagonal().real();

                //eigen_sparse_real_t overlap_root = (mat_current->basis().dot(mat_previous->basis())).pruned(1e-4,0.5).cwiseAbs();
                //eigen_sparse_real_t overlap_root = (mat_previous->basis()).cwiseAbs(); // hiermit waere es viel schneller

                /*eigen_sparse_t overlap_tmp(mat_current->basis().innerSize(),mat_current->basis().innerSize());
                std::vector<eigen_triplet_t> triplets;
                triplets.reserve(std::max(mat_current->basis().nonZeros(),mat_previous->basis().nonZeros())*mat_current->basis().outerSize());
                //std::cout << 1000 << std::endl;
                //scalar_t number = 0;
                for (eigen_idx_t k=0; k<mat_current->basis().outerSize(); ++k) {
                    for (eigen_iterator_t it_i(mat_current->basis(),k); it_i; ++it_i) {
                        for (eigen_iterator_t it_j(mat_previous->basis(),k); it_j; ++it_j) {
                            //number += it_i.value()*it_j.value();
                            triplets.push_back(eigen_triplet_t(it_i.index(),it_j.index(),0));//std::conj(it_i.value())*it_j.value()));
                        }
                    }
                }
                //std::cout << 1001 << std::endl;
                overlap_tmp.setFromTriplets(triplets.begin(), triplets.end());
                //std::cout << 1002 << std::endl;
                eigen_sparse_real_t overlap_root = (mat_previous->basis()).cwiseAbs();
                //std::cout << 1003 << std::endl;*/

                std::unordered_map<eigen_idx_t,eigen_idx_t> lines_idx_new;

                for (eigen_idx_t k=0; k<overlap_root.outerSize(); ++k) {
                    eigen_idx_t idx_current = k;

                    eigen_idx_t idx_previous = -1;
                    for (eigen_iterator_real_t it(overlap_root,k); it; ++it) {
                        if (it.value() > std::sqrt(0.5)) {
                            idx_previous = it.index();
                            break;
                        }
                    }

                    real_t eigenenergies_current = eigenenergies[idx_current];
                    if (idx_previous >= 0) {
                        //lines[lines_idx[idx_previous]].push_back(eigenenergies_current);
                        triplets_lines_eigenenergies.push_back(eigen_triplet_real_t(lines_idx[idx_previous],step,eigenenergies_current));
                        lines_idx_new[idx_current] = lines_idx[idx_previous];
                    } else {
                        //lines.push_back(std::vector<real_t>(1,eigenenergies_current));
                        triplets_lines_eigenenergies.push_back(eigen_triplet_real_t(lines_idx_max,step,eigenenergies_current));
                        lines_idx_new[idx_current] = lines_idx_max++;
                    }
                }

                lines_idx = lines_idx_new;
                mat_previous = mat_current;
                step++;
            }

            eigen_sparse_real_t lines_eigenenergies(lines_idx_max,step);
            lines_eigenenergies.setFromTriplets(triplets_lines_eigenenergies.begin(), triplets_lines_eigenenergies.end());

            // save lines
            byte_t lines_flags = 0;
            if (lines_eigenenergies.IsRowMajor) {
                lines_flags |= csr_not_csc;
            }
            storage_idx_t lines_rows = lines_eigenenergies.rows();
            storage_idx_t lines_cols = lines_eigenenergies.cols();
            std::vector<storage_real_t> lines_data(lines_eigenenergies.valuePtr(), lines_eigenenergies.valuePtr()+lines_eigenenergies.nonZeros());
            std::vector<storage_idx_t> lines_indices(lines_eigenenergies.innerIndexPtr(), lines_eigenenergies.innerIndexPtr()+lines_eigenenergies.nonZeros());
            std::vector<storage_idx_t> lines_indptr(lines_eigenenergies.outerIndexPtr(), lines_eigenenergies.outerIndexPtr()+lines_eigenenergies.outerSize());

            bytes_t bytes;
            Serializer s;
            s << lines_flags;
            s << lines_rows;
            s << lines_cols;
            s << lines_data;
            s << lines_indices;
            s << lines_indptr;
            s.save(bytes);

            // open file
            FILE *pFile;
            pFile = fopen("output/lines.mat" , "wb" ); // filename_.c_str() // TODO

            // write
            fwrite(&bytes[0], 1 , sizeof(byte_t)*bytes.size(), pFile );

            // close file
            fclose(pFile);
        }
    }

protected:
    std::shared_ptr<Hamiltonianmatrix> doProcessing(std::shared_ptr<Hamiltonianmatrix> work) {
        // if results can not be loaded
        if (!work->exist() || !work->load()) {

            if (work->num_basisvectors() > 1) {
                // diagonalization
                Eigen::SelfAdjointEigenSolver<eigen_dense_t> eigensolver(eigen_dense_t(work->entries()));

                // eigenvalues and eigenvectors
                eigen_vector_real_t evals = eigensolver.eigenvalues();
                eigen_sparse_t evecs = eigensolver.eigenvectors().sparseView(1e-4,0.5);

                work->entries().setZero();
                work->entries().reserve(evals.size());
                for (eigen_idx_t idx = 0; idx < evals.size(); ++idx) {
                    work->entries().insert(idx, idx) = evals.coeffRef(idx);
                }
                work->entries().makeCompressed();

                work->basis() = (work->basis() * evecs).pruned(1e-4,0.5);
            }

            // save result
            work->save();
        }

        return work;
    }

    void doPrework(size_t numWork) {
        std::cout << ">>DIM" << std::setw(7) << matrix_dimension[numWork] << std::endl;
    }

    void doPostwork(size_t numWork) {
        std::cout << ">>OUT" << std::setw(7) << numWork+1 << std::setw(7) << matrix_step[numWork] << std::setw(7) << matrix_blocks[numWork] << std::setw(7) << matrix_block[numWork] << " " << matrix_path[numWork] << std::endl;
    }

    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix;
    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix_diag; // TODO figure out whether the same pointers as in matrix
    std::vector<std::shared_ptr<Configuration>> params;
    std::vector<std::string> matrix_path;
    std::vector<size_t> matrix_step;
    std::vector<size_t> matrix_blocks;
    std::vector<size_t> matrix_block;
    std::vector<size_t> matrix_dimension;
};

class HamiltonianOne : public Hamiltonian{
public:
    HamiltonianOne(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<MpiEnvironment> mpi, std::shared_ptr<BasisnamesOne> basis_one) : Hamiltonian(mpi), basis_one(basis_one), path_cache(path_cache) {
        configure(config);
        build();
    }

    std::shared_ptr<const BasisnamesOne> names() const {
        return basis_one;
    }

    const Configuration& getConf() const { // TODO in Configurable Klasse auslagern, von der geerbt werrden soll
        return conf;
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
        conf = basis_one->getConf();
        conf["deltaESingle"] = config["deltaESingle"];
        conf["diamagnetism"] = config["diamagnetism"];

        conf["deltaESingle"] >> deltaE;
        conf["species1"] >> species;

        diamagnetism = conf["diamagnetism"].str() == "true";

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

        if (mpi->rank() == 0) {
            real_t tol = 1e-32;

            /*real_t energycutoff = 0.7e-5;

            real_t theta = 0*M_PI/180.;
            real_t theta2 = 60*M_PI/180.;

            real_t min_E_x = 0;
            real_t min_E_y = 0;
            real_t min_E_z = 0;
            real_t max_E_x = 0;
            real_t max_E_y = 1e-11*std::sin(theta2);
            real_t max_E_z = 1e-11*std::cos(theta2);
            real_t min_B_x = 0;
            real_t min_B_y = 50*4.254382e-10*std::sin(theta);
            real_t min_B_z = 50*4.254382e-10*std::cos(theta);
            real_t max_B_x = 0;
            real_t max_B_y = 50*4.254382e-10*std::sin(theta);
            real_t max_B_z = 50*4.254382e-10*std::cos(theta);*/

            // === calculate one-atom Hamiltonians ===
            scalar_t min_E_0, min_E_p, min_E_m, min_B_0, min_B_p, min_B_m, max_E_0, max_E_p, max_E_m, max_B_0, max_B_p, max_B_m;
            changeToSpherical(min_E_x, min_E_y, min_E_z, min_E_p, min_E_m, min_E_0);
            changeToSpherical(max_E_x, max_E_y, max_E_z, max_E_p, max_E_m, max_E_0);
            changeToSpherical(min_B_x, min_B_y, min_B_z, min_B_p, min_B_m, min_B_0);
            changeToSpherical(max_B_x, max_B_y, max_B_z, max_B_p, max_B_m, max_B_0);

            bool exist_E_0 = (std::abs(min_E_0) != 0 || std::abs(max_E_0) != 0);
            bool exist_E_1 = (std::abs(min_E_p) != 0 || std::abs(max_E_p) != 0);
            bool exist_B_0 = (std::abs(min_B_0) != 0 || std::abs(max_B_0) != 0);
            bool exist_B_1 = (std::abs(min_B_p) != 0 || std::abs(max_B_p) != 0);

            // --- count entries of Hamiltonian parts ---
            size_t size_basis = basis_one->size();
            size_t size_energy = basis_one->size();

            // --- construct energy Hamiltonian part ---
            std::cout << "One-atom Hamiltonian, construct diagonal hamiltonian" << std::endl;

            Hamiltonianmatrix hamiltonian_energy(size_basis, size_energy);

            real_t energy_initial = 0;
            for (const auto &state: basis_one->initial()) {
                energy_initial += energy_level(species,state.n,state.l,state.j);
            }
            energy_initial /= basis_one->initial().size();

            std::vector<bool> is_necessary(basis_one->size(),false);
            idx_t idx = 0;
            for (const auto &state : *basis_one) {
                real_t val = energy_level(species,state.n,state.l,state.j)-energy_initial;
                if (std::abs(val) < deltaE+1e-11 || deltaE < 0) { // TODO
                    is_necessary[state.idx] = true;
                    hamiltonian_energy.addEntries(idx,idx,val);
                    hamiltonian_energy.addBasis(idx,idx,1);
                    ++idx;
                }
            }
            std::cout << "One-atom basis, size without restrictions in energy: " << basis_one->size() << std::endl;

            basis_one->removeUnnecessaryStates(is_necessary);

            hamiltonian_energy.compress(basis_one->dim(), basis_one->dim());

            std::cout << "One-atom basis, size with restrictions in energy: " << basis_one->size() << std::endl;



            std::cout << ">>BAS" << std::setw(7) << basis_one->size() << std::endl;

            // initialize uuid generator
            boost::uuids::random_generator generator;

            // generate uuid
            std::string uuid;
            boost::uuids::uuid u = generator();
            boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

            // save basis
            boost::filesystem::path path_basis = boost::filesystem::temp_directory_path();
            path_basis /= "basis_one_"+uuid+".csv";
            basis_one->save(path_basis.string());

            std::cout << ">>STA " << path_basis.string() << std::endl;

            // --- precalculate matrix elements ---
            std::cout << "One-atom Hamiltonian, precalculate matrix elements" << std::endl;

            MatrixElements matrix_elements(conf, species, (path_cache / "cache_elements.db").string());

            if (exist_E_0) matrix_elements.precalculateElectricMomentum(basis_one, 0);
            if (exist_E_1) matrix_elements.precalculateElectricMomentum(basis_one, 1);
            if (exist_E_1) matrix_elements.precalculateElectricMomentum(basis_one, -1);

            if (exist_B_0) matrix_elements.precalculateMagneticMomentum(basis_one, 0);
            if (exist_B_1) matrix_elements.precalculateMagneticMomentum(basis_one, 1);
            if (exist_B_1) matrix_elements.precalculateMagneticMomentum(basis_one, -1);

            if (diamagnetism && (exist_B_0 || exist_B_1)) matrix_elements.precalculateDiamagnetism(basis_one, 0, 0);
            if (diamagnetism && (exist_B_0 || exist_B_1)) matrix_elements.precalculateDiamagnetism(basis_one, 2, 0);
            if (diamagnetism && exist_B_0 && exist_B_1) matrix_elements.precalculateDiamagnetism(basis_one, 2, 1);
            if (diamagnetism && exist_B_0 && exist_B_1) matrix_elements.precalculateDiamagnetism(basis_one, 2, -1);
            if (diamagnetism && exist_B_1) matrix_elements.precalculateDiamagnetism(basis_one, 2, 2);
            if (diamagnetism && exist_B_1) matrix_elements.precalculateDiamagnetism(basis_one, 2, -2);

            // --- count entries of Hamiltonian parts ---
            std::cout << "One-atom Hamiltonian, count number of entries within the field Hamiltonian" << std::endl;

            size_basis = basis_one->size();
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

            for (const auto &state_col : *basis_one) {
                for (const auto &state_row : *basis_one) {
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

            // --- construct field Hamiltonian parts ---
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

            for (const auto &state_col : *basis_one) {
                for (const auto &state_row : *basis_one) {
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

            hamiltonian_electricMomentum_0.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_electricMomentum_p.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_electricMomentum_m.compress(basis_one->dim(), basis_one->dim());

            hamiltonian_magneticMomentum_0.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_magneticMomentum_p.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_magneticMomentum_m.compress(basis_one->dim(), basis_one->dim());

            hamiltonian_diamagnetism_00.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_diamagnetism_20.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_diamagnetism_2p.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_diamagnetism_2m.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_diamagnetism_2pp.compress(basis_one->dim(), basis_one->dim());
            hamiltonian_diamagnetism_2mm.compress(basis_one->dim(), basis_one->dim());

            // --- construct Hamiltonians ---
            std::cout << "One-atom Hamiltonian, assemble Hamiltonians" << std::endl;

            matrix.reserve(nSteps);
            params.reserve(nSteps);

            // open database
            boost::filesystem::path path_db;

            if (utils::is_complex<scalar_t>::value) {
                path_db = path_cache / "cache_matrix_complex.db";
            } else {
                path_db = path_cache / "cache_matrix_real.db";
            }
            SQLite3 db(path_db.string());

            // loop through steps
            for (size_t step = 0; step < nSteps; ++step) {
                real_t normalized_position = (nSteps > 1) ? step/(nSteps-1.) : 0;

                std::shared_ptr<Configuration> par = std::make_shared<Configuration>(conf);

                // save fields to ConfParser object
                (*par)["Ex"] = min_E_x+normalized_position*(max_E_x-min_E_x);
                (*par)["Ey"] = min_E_y+normalized_position*(max_E_y-min_E_y);
                (*par)["Ez"] = min_E_z+normalized_position*(max_E_z-min_E_z);
                (*par)["Bx"] = min_B_x+normalized_position*(max_B_x-min_B_x);
                (*par)["By"] = min_B_y+normalized_position*(max_B_y-min_B_y);
                (*par)["Bz"] = min_B_z+normalized_position*(max_B_z-min_B_z);

                // create table if necessary
                std::stringstream query;
                std::string spacer = "";

                if (step == 0) {
                    query << "CREATE TABLE IF NOT EXISTS cache_one (uuid text NOT NULL PRIMARY KEY, "
                             "created TIMESTAMP DEFAULT CURRENT_TIMESTAMP, "
                             "accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP";
                    for (auto p: *par) {
                        query << ", " << p.key << " text";
                    }
                    query << ", UNIQUE (";
                    for (auto p: *par) {
                        query << spacer << p.key;
                        spacer = ", ";
                    }
                    query << "));";
                    db.exec(query.str());
                }

                // get uuid as filename
                std::string uuid;

                query.str(std::string());
                spacer = "";
                query << "SELECT uuid FROM cache_one WHERE ";
                for (auto p: *par) {
                    query << spacer << p.key << "='" << p.value.str() << "'";
                    spacer = " AND ";
                }
                query << ";";
                SQLite3Result result = db.query(query.str());

                if (result.size() == 1) {
                    uuid = result.first()->str();

                    /*query.str(std::string());
                    query << "UPDATE cache_one SET accessed = CURRENT_TIMESTAMP WHERE uuid = '" << uuid << "';";
                    db.exec(query.str());*/ // TODO This is very slow on rqo-donkey!
                } else {
                    boost::uuids::uuid u = generator();
                    boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

                    query.str(std::string());
                    query << "INSERT INTO cache_one (uuid";
                    for (auto p: *par) {
                        query << ", " << p.key;
                    }
                    query << ") values ( '" << uuid << "'";
                    for (auto p: *par) {
                        query << ", " << "'" << p.value.str() << "'";
                    }
                    query << ");";
                    db.exec(query.str());
                }

                // check whether .mat and .json file exists and compare settings in program with settings in .json file
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
                        if (*par == params_loaded) {
                            is_existing = true;
                        }
                    }
                }

                // create .json file if "is_existing" is false
                if (!is_existing) {
                    par->save_to_json(path_json.string());
                }

                // calculate Hamiltonian if "is_existing" is false
                std::shared_ptr<Hamiltonianmatrix> mat;
                if (!is_existing) {
                    scalar_t E_0 = min_E_0+normalized_position*(max_E_0-min_E_0);
                    scalar_t E_p = min_E_p+normalized_position*(max_E_p-min_E_p);
                    scalar_t E_m = min_E_m+normalized_position*(max_E_m-min_E_m);
                    scalar_t B_0 = min_B_0+normalized_position*(max_B_0-min_B_0);
                    scalar_t B_p = min_B_p+normalized_position*(max_B_p-min_B_p);
                    scalar_t B_m = min_B_m+normalized_position*(max_B_m-min_B_m);

                    mat = std::make_shared<Hamiltonianmatrix>(hamiltonian_energy
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
                                                              -std::sqrt(1.5)*hamiltonian_diamagnetism_2mm*B_p*B_p
                                                              );
                } else {
                    mat = std::make_shared<Hamiltonianmatrix>();
                    //mat->compress(hamiltonian_energy.num_basisvectors(),hamiltonian_energy.num_coordinates());
                }

                // save everything
                mat->addFilename(path_mat.string());
                mat->addIsExisting(is_existing);
                matrix.push_back(std::move(mat));
                params.push_back(std::move(par));
                matrix_path.push_back(path.string());
                matrix_step.push_back(step);
                matrix_blocks.push_back(1); // TODO
                matrix_block.push_back(0);
                matrix_dimension.push_back(hamiltonian_energy.num_basisvectors());

                std::cout << "One-atom Hamiltonian, " <<  matrix.size() << ". Hamiltonian assembled" << std::endl;
            }

            std::cout << ">>TOT" << std::setw(7) << matrix.size() << std::endl;

            std::cout << "One-atom Hamiltonian, all Hamiltonians assembled" << std::endl;
        }

        // === diagonalize matrices using MpiLoadbalancingSimple ===
        run(matrix, matrix_diag);
    }

private:
    Configuration conf;
    std::shared_ptr<BasisnamesOne> basis_one;
    real_t deltaE;
    real_t min_E_x,min_E_y,min_E_z,max_E_x,max_E_y,max_E_z,min_B_x,min_B_y,min_B_z,max_B_x,max_B_y,max_B_z;
    size_t nSteps;
    bool diamagnetism;
    std::string species;
    boost::filesystem::path path_cache;

};

class HamiltonianTwo : public Hamiltonian {
public:
    HamiltonianTwo(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<MpiEnvironment> mpi, std::shared_ptr<const HamiltonianOne> hamiltonian_one)  :
        Hamiltonian(mpi), hamiltonian_one1(hamiltonian_one), hamiltonian_one2(hamiltonian_one), basis_two(std::make_shared<BasisnamesTwo>(hamiltonian_one->names())), path_cache(path_cache) { // TODO

        samebasis = true;


        /*size_t nSteps_one = hamiltonian_one1->size(); // TODO

        // TODO : in extra methode auslagern, hamiltonian_one2 mitberuecksichtigen (alle delta in externes Objekt auslagern)
        conf = basis_two->getConf();
        conf["deltaE"] = config["deltaE"];

        conf["deltaE"] >> deltaE;
        deltaE = std::fmax(deltaE,1e-24); // TODO remove hack
        conf["species1"] >> species1;
        conf["species2"] >> species2;

        config["steps"] >> nSteps_two;
        config["minR"] >> min_R;
        config["maxR"] >> max_R;

        dipoledipole = config["dd"].str() == "true"; // TODO

        if (min_R == max_R && nSteps_one == 1){
            nSteps_two = 1;
        } else {
            config["steps"] >> nSteps_two;
        }*/

        calculate(config);
    }

    HamiltonianTwo(const Configuration &config, boost::filesystem::path& path_cache, std::shared_ptr<MpiEnvironment> mpi, std::shared_ptr<const HamiltonianOne> hamiltonian_one1, std::shared_ptr<const HamiltonianOne> hamiltonian_one2) :
        Hamiltonian(mpi), hamiltonian_one1(hamiltonian_one1), hamiltonian_one2(hamiltonian_one2), basis_two(std::make_shared<BasisnamesTwo>(hamiltonian_one1->names(), hamiltonian_one2->names())), path_cache(path_cache) {

        samebasis = false;

        /*size_t nSteps_one = hamiltonian_one1->size(); // TODO

        // TODO : in extra methode auslagern, hamiltonian_one2 mitberuecksichtigen (alle delta in externes Objekt auslagern)
        conf = basis_two->getConf();
        conf["deltaE"] = config["deltaE"];

        conf["deltaE"] >> deltaE;
        deltaE = std::fmax(deltaE,1e-24); // TODO remove hack
        conf["species1"] >> species1;
        conf["species2"] >> species2;

        config["steps"] >> nSteps_two;
        config["minR"] >> min_R;
        config["maxR"] >> max_R;

        dipoledipole = config["dd"].str() == "true"; // TODO

        if (min_R == max_R && nSteps_one == 1){
            nSteps_two = 1;
        } else {
            config["steps"] >> nSteps_two;
        }*/

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

        if (mpi->rank() == 0) {
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
            Configuration conf_matpair = basis_two->getConf();
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
                //conf_mat.push_back(conf_matsingle + conf_matpair); // conf_matpair overwrites settings in conf_matsingle // TODO !!!!!!!!!!!!!!!!!!
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



            // --- determine coordinates whose usage is not restricted ---
            std::cout << "Two-atom basis, size without restrictions: " << basis_two->size() << std::endl;
            std::cout << "Two-atom basis, determine pair states whose usage is not forbidden because of restrictions in energy or quatum numbers" << std::endl;

            auto basis_one1 = hamiltonian_one1->names();
            auto basis_one2 = hamiltonian_one2->names();
            std::vector<bool> notrestricted_coordinate1(basis_one1->size(), false);
            std::vector<bool> notrestricted_coordinate2(basis_one2->size(), false);
            std::vector<bool> notrestricted_coordinate(basis_two->size(), false);
            std::vector<StateOne> initial1 = basis_one1->initial();
            std::vector<StateOne> initial2 = basis_one2->initial();
            StateTwo initial = basis_two->initial();

            // coordinates of atom 1
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
                    notrestricted_coordinate1[state.idx] = true;
                }
            }

            // coordinates of atom 2
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
                    notrestricted_coordinate2[state.idx] = true;
                }
            }

            // coordinates of the atom pair
            float M = initial.m[0]+initial.m[1];
            int parity = (initial.l[0]+initial.l[1]) % 2;

            for (const auto &state: *basis_two) {
                if (conserveM && state.m[0]+state.m[1] != M) continue;
                if (conserveParityL && (state.l[0]+state.l[1]) % 2 != parity) continue;

                notrestricted_coordinate[state.idx] = true;
            }


            // === construct pair Hamiltonian consistent of combined one-atom Hamiltonians (1 x Hamiltonian2 + Hamiltonian1 x 1)  ===


            // --- determine necessary symmetries ---
            std::cout << "Two-atom basis, determine symmetrized subspaces" << std::endl;
            enum symmetries_t {ALL, SYM, ASYM};

            std::vector<symmetries_t> symmetries;
            if (samebasis && multipoleexponent <= 3) {
                symmetries.push_back(SYM);
                if (initial.first() != initial.second()) {
                    symmetries.push_back(ASYM);
                }
            } else {
                symmetries.push_back(ALL);
            }


            // --- load one-atom Hamiltonians ---
            std::cout << "Two-atom Hamiltonian, construct contribution of combined one-atom Hamiltonians" << std::endl;

            std::map<symmetries_t,std::vector<Hamiltonianmatrix>> mat_single;
            for (symmetries_t symmetry : symmetries) {
                mat_single[symmetry].reserve(nSteps_one);
            }

            // combine the Hamiltonians of the two atoms, beeing aware of the energy cutoff and the list of necessary coordinates
            for (size_t i = 0; i < nSteps_one; ++i) {
                for (symmetries_t symmetry : symmetries) {
                    if (symmetry == SYM) mat_single[symmetry].push_back(combineSym(*(hamiltonian_one1->get(i)), *(hamiltonian_one2->get(i)), deltaE, notrestricted_coordinate1, notrestricted_coordinate2, notrestricted_coordinate));
                    if (symmetry == ASYM) mat_single[symmetry].push_back(combineAsym(*(hamiltonian_one1->get(i)), *(hamiltonian_one2->get(i)), deltaE, notrestricted_coordinate1, notrestricted_coordinate2, notrestricted_coordinate));
                    if (symmetry == ALL) mat_single[symmetry].push_back(combineAll(*(hamiltonian_one1->get(i)), *(hamiltonian_one2->get(i)), deltaE, notrestricted_coordinate1, notrestricted_coordinate2, notrestricted_coordinate));
                }
                std::cout << "Two-atom Hamiltonian, "<< i+1 << ". Hamiltonian combined" << std::endl;
            }

            // --- remove unnecessary (i.e. more or less empty) basisvectors ---
            std::cout << "Two-atom basis, remove unnecessary basis vectors (i.e. vectors with too small norm)" << std::endl;

            for (size_t i = 0; i < nSteps_one; ++i) {
                for (symmetries_t symmetry : symmetries) {
                    mat_single[symmetry][i].removeUnnecessaryBasisvectors();
                }
            }

            // --- find coordinates that are used ---
            std::cout << "Two-atom basis, remove unnecessary states (i.e. states with too unlikely occurence)" << std::endl;
            std::vector<bool> used_coordinate(basis_two->size(), false);

            for (size_t i = 0; i < nSteps_one; ++i) {
                for (symmetries_t symmetry : symmetries) {
                    mat_single[symmetry][i].findUnnecessaryStates(used_coordinate);
                }
            }

            // --- remove coordinates that are not used --- // TODO make this not necessary by restricting the coordinates to sensible ones starting from the beginning (this remark is connected to all other remarks marked with [*])
            for (size_t i = 0; i < nSteps_one; ++i) {
                for (symmetries_t symmetry : symmetries) {
                    mat_single[symmetry][i].removeUnnecessaryStates(used_coordinate);
                }
            }
            basis_two->removeUnnecessaryStates(used_coordinate);


            /*std::cout << mat_single_sym[0].num_coordinates() << std::endl; // TODO
            std::cout << mat_single_sym.back().num_basisvectors() << std::endl;// TODO
            std::cout << mat_single_sym[0].num_basisvectors() << std::endl;// TODO*/

            std::cout << "Two-atom basis, number of necessary states: " << basis_two->size() << std::endl;

            std::cout << ">>BAS" << std::setw(7) << basis_two->size() << std::endl;

            // --- find index of initial state ---

            // initialize uuid generator
            boost::uuids::random_generator generator;

            // generate uuid
            std::string uuid;
            boost::uuids::uuid u = generator();
            boost::algorithm::hex(u.begin(), u.end(), std::back_inserter(uuid));

            // save basis
            boost::filesystem::path path_basis = boost::filesystem::temp_directory_path();
            path_basis /= "basis_two_"+uuid+".csv";
            basis_two->save(path_basis.string());

            std::cout << ">>STA " << path_basis.string() << std::endl;

            // === construct pair Hamiltonians for all orders of the multipole expansion  ===

            std::vector<int> exponent_multipole;
            std::vector<Hamiltonianmatrix> mat_multipole;
            MatrixElements matrixelements_atom1(conf_tot, species1, (path_cache / "cache_elements.db").string());
            MatrixElements matrixelements_atom2(conf_tot, species2, (path_cache / "cache_elements.db").string());
            std::vector<idx_t> size_mat_multipole;

            int idx_multipole_max = -1;

            if (multipoleexponent > 2) {

                // --- initialize two-atom interaction Hamiltonians ---
                std::cout << "Two-atom hamiltonian, initialize interaction Hamiltonians" << std::endl;

                int kappa_min = 1; // spherical dipole operators
                int kappa_max = multipoleexponent-kappa_min-1;
                int sumOfKappas_min = kappa_min+kappa_min;
                int sumOfKappas_max = kappa_max+kappa_min;
                idx_multipole_max = sumOfKappas_max-sumOfKappas_min;

                exponent_multipole.reserve(idx_multipole_max+1);
                mat_multipole.reserve(idx_multipole_max+1);
                size_mat_multipole.resize(idx_multipole_max+1);

                // --- precalculate matrix elements ---
                std::cout << "Two-atom basis, get one-atom states needed for the two-atom basis"<< std::endl;

                basis_one1 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromFirst(basis_two));
                basis_one2 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromSecond(basis_two));

                for (int kappa = kappa_min; kappa<=kappa_max; ++kappa) {
                    std::cout << "Two-atom hamiltonian, precalculate matrix elements for kappa = " << kappa << std::endl;
                    matrixelements_atom1.precalculateMultipole(basis_one1, kappa);
                    matrixelements_atom2.precalculateMultipole(basis_one2, kappa);
                }

                // --- count entries of two-atom interaction Hamiltonians ---
                std::cout << "Two-atom Hamiltonian, count number of entries within the interaction Hamiltonians" << std::endl;

                for (int sumOfKappas = sumOfKappas_min; sumOfKappas<=sumOfKappas_max; ++sumOfKappas) {
                    int idx_multipole = sumOfKappas-sumOfKappas_min;

                    for (const auto &state_col : *basis_two) {
                        //if (!used_coordinate[state_col.idx]) continue; // TODO (this remark is connected to all other remarks marked with [*])
                        int M_col = state_col.first().m + state_col.second().m;

                        for (const auto &state_row : *basis_two) {
                            //if (!used_coordinate[state_row.idx]) continue; // TODO (this remark is connected to all other remarks marked with [*])
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

                // --- construct two-atom interaction Hamiltonians ---
                size_t size_basis = basis_two->size();

                for (int sumOfKappas = sumOfKappas_min; sumOfKappas<=sumOfKappas_max; ++sumOfKappas) {
                    std::cout << "Two-atom Hamiltonian, construct interaction Hamiltonian that belongs to 1/R^" << sumOfKappas+1 << std::endl;

                    int idx_multipole = sumOfKappas-sumOfKappas_min;

                    exponent_multipole.push_back(sumOfKappas+1);
                    mat_multipole.push_back(Hamiltonianmatrix(size_basis, 2*size_mat_multipole[idx_multipole])); // factor of 2 because triangular matrix is not sufficient

                    for (const auto &state_col : *basis_two) {
                        //if (!used_coordinate[state_col.idx]) continue; // TODO (this remark is connected to all other remarks marked with [*])
                        int M_col = state_col.first().m + state_col.second().m;

                        for (const auto &state_row : *basis_two) {
                            //if (!used_coordinate[state_row.idx]) continue; // TODO (this remark is connected to all other remarks marked with [*])
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

                    mat_multipole[idx_multipole].compress(basis_two->dim(), basis_two->dim()); // TODO substitute dim() by size()
                }
            }

            // === construct Hamiltonians === // TODO Logik in eigene Klasse
            std::cout << "Two-atom Hamiltonian, assemble Hamiltonians" << std::endl;

            matrix.reserve(nSteps_two);

            // open database
            boost::filesystem::path path_db;

            if (utils::is_complex<scalar_t>::value) {
                path_db = path_cache / "cache_matrix_complex.db";
            } else {
                path_db = path_cache / "cache_matrix_real.db";
            }
            SQLite3 db(path_db.string());

            // initialize variables
            Configuration conf;
            size_t step_one = 0;

            std::map<symmetries_t,std::vector<Hamiltonianmatrix>> mat_multipole_transformed;
            for (symmetries_t symmetry : symmetries) {
                mat_multipole_transformed[symmetry].resize(idx_multipole_max+1);
            }

            bool flag_perhapsmissingtable = true;

            std::map<symmetries_t,std::string> symmetries_name;
            symmetries_name[SYM] = "sym";
            symmetries_name[ASYM] = "asym";
            symmetries_name[ALL] = "all";

            std::vector<Hamiltonianmatrix> blocks_mat_single;
            std::vector<std::vector<Hamiltonianmatrix>> blocks_mat_multipole;
            std::vector<std::string> blocks_symmetries_name;
            std::vector<size_t> blocks_identification;
            blocks_mat_single.reserve(symmetries.size());
            blocks_mat_multipole.reserve(symmetries.size());
            blocks_symmetries_name.reserve(symmetries.size());
            blocks_identification.reserve(symmetries.size());

            // TODO make this not necessary (this remark is connected to all other remarks marked with [*])
            std::string symmetryies_string = "";
            for (symmetries_t symmetry : symmetries) {
                if (symmetryies_string != "") symmetryies_string += " ";
                symmetryies_string += symmetries_name[symmetry];
            }

            // loop through steps
            for (size_t step_two = 0; step_two < nSteps_two; ++step_two) {
                real_t normalized_position = (nSteps_two > 1) ? step_two/(nSteps_two-1.) : 0;
                real_t position = min_R+normalized_position*(max_R-min_R);

                if (step_two == 0 || (nSteps_one > 1 && step_one <= step_two)) {
                    conf = conf_mat[step_one];

                    blocks_mat_multipole.clear();
                    blocks_mat_single.clear();
                    blocks_symmetries_name.clear();
                    blocks_identification.clear();

                    for (symmetries_t symmetry : symmetries) {
                        Hamiltonianmatrix abstotalmatrix = mat_single[symmetry][step_one];
                        std::vector<Hamiltonianmatrix> transformedmatrix;
                        transformedmatrix.reserve(idx_multipole_max+1);
                        for (int idx_multipole = 0; idx_multipole <= idx_multipole_max; ++idx_multipole) {
                            transformedmatrix.push_back(mat_multipole[idx_multipole].changeBasis(mat_single[symmetry][step_one].basis()));
                            abstotalmatrix += transformedmatrix.back().abs();
                        }

                        // TODO find submatrices
                        StateTwo initial_state = basis_two->initial();
                        std::set<ptrdiff_t> initial_indices;
                        if (conserveM) {
                            initial_indices.insert(initial_state.idx);
                        } else {
                            for (auto state : *basis_two) {
                                if ((state.n[0] == initial_state.n[0]) && (state.l[0] == initial_state.l[0])  && (state.j[0] == initial_state.j[0])  &&
                                        (state.n[1] == initial_state.n[1]) && (state.l[1] == initial_state.l[1])  && (state.j[1] == initial_state.j[1])) {
                                    initial_indices.insert(state.idx);
                                }
                            }
                        }

                        std::vector<bool> isNecessaryBasisvector(abstotalmatrix.basis().cols(),false);
                        for (eigen_idx_t k_1=0; k_1<abstotalmatrix.basis().outerSize(); ++k_1) {
                            for (eigen_iterator_t triple(abstotalmatrix.basis(),k_1); triple; ++triple) {
                                ptrdiff_t col = triple.col(); // basis vector
                                ptrdiff_t row = triple.row(); // coordinate
                                if (initial_indices.count(row)) {
                                    isNecessaryBasisvector[col] = true;
                                }
                            }
                        }

                        std::vector<size_t> idxNecessaryBasisvector;
                        for (size_t i = 0; i < isNecessaryBasisvector.size(); ++i) {
                            if (isNecessaryBasisvector[i]) {
                                idxNecessaryBasisvector.push_back(i);
                            }
                        }



                        /*eigen_sparse_t tmp = abstotalmatrix.entries();
                        for (int i = 0; i < 4; ++i) {
                            tmp = tmp*tmp;
                        }

                        std::vector<std::vector<ptrdiff_t>> blockIndices;
                        for (auto& i: idxNecessaryBasisvector) {
                            // TODO combine if one of the indices agree
                            std::cout << bas.cols() << " " << i << std::endl;
                            blockIndices.push_back(std::vector<ptrdiff_t>());
                            for (eigen_iterator_t triple(tmp,i); triple; ++triple) {
                                std::cout << triple.index() << " " << triple.value() << std::endl;
                                blockIndices.back().push_back(triple.index());
                            }
                        }*/

                        std::vector<std::set<size_t>> blockSets; // TODO do not check for blocks if no interaction
                        for (eigen_idx_t k_1=0; k_1<abstotalmatrix.entries().outerSize(); ++k_1) {

                            std::vector<size_t> current_indices;
                            for (eigen_iterator_t triple(abstotalmatrix.entries(),k_1); triple; ++triple) {
                                current_indices.push_back(triple.index());
                            }

                            bool breakloop = false;
                            for (auto& set: blockSets) { // TODO make function out of it and use return
                                for (auto& i: current_indices) {
                                    if (set.count(i)) {
                                        breakloop = true;
                                        break;
                                    }
                                }
                                if (breakloop) {
                                    set.insert(current_indices.begin(), current_indices.end());
                                    break;
                                }
                            }
                            if (!breakloop) {
                                blockSets.push_back(std::set<size_t>(current_indices.begin(), current_indices.end()));
                            }
                        }

                        std::vector<std::vector<ptrdiff_t>> blockIndices;
                        blockIndices.reserve(blockSets.size());
                        bool nottrivial = false;
                        for (auto& set: blockSets) {
                            for (auto& i: idxNecessaryBasisvector) {
                                if (set.count(i)) {
                                    blockIndices.push_back(std::vector<ptrdiff_t>(set.begin(), set.end()));
                                    nottrivial |= set.size() > 1;
                                    break;
                                }
                            }
                        }
                        if (!nottrivial) {
                            blockIndices.clear();
                            blockIndices.push_back(std::vector<ptrdiff_t>(abstotalmatrix.num_basisvectors()));
                            std::iota (blockIndices.back().begin(), blockIndices.back().end(), 0);
                        }

                        for (auto& indices: blockIndices) {
                            blocks_mat_multipole.push_back(std::vector<Hamiltonianmatrix>());
                            blocks_mat_multipole.back().reserve(idx_multipole_max+1);
                            for (int idx_multipole = 0; idx_multipole <= idx_multipole_max; ++idx_multipole) {
                                blocks_mat_multipole.back().push_back(transformedmatrix[idx_multipole].getBlock(indices));
                            }
                            blocks_mat_single.push_back(mat_single[symmetry][step_one].getBlock(indices));
                            blocks_symmetries_name.push_back(symmetries_name[symmetry]);

                            size_t identification = 0;
                            for (auto& i: indices) utils::hash_combine(identification, i);
                            blocks_identification.push_back(identification);
                        }
                    }

                    ++step_one;
                }

                conf["R"] = position;
                conf["allsymmetries"] = symmetryies_string; // TODO make this not necessary (this remark is connected to all other remarks marked with [*])

                // loop through symmetries
                for (size_t idx_block = 0; idx_block < blocks_mat_single.size(); ++idx_block) {
                    conf["symmetry"] = blocks_symmetries_name[idx_block];
                    conf["sub"] = blocks_identification[idx_block];

                    Hamiltonianmatrix totalmatrix = blocks_mat_single[idx_block];
                    for (int idx_multipole = 0; idx_multipole <= idx_multipole_max; ++idx_multipole) {
                        real_t pos = 1./std::pow(position,exponent_multipole[idx_multipole]);
                        totalmatrix += blocks_mat_multipole[idx_block][idx_multipole]*pos;
                    }

                    // create table if necessary
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
                        db.exec(query.str());

                        flag_perhapsmissingtable = false;
                    }

                    // get uuid as filename
                    std::string uuid;

                    query.str(std::string());
                    spacer = "";
                    query << "SELECT uuid FROM cache_two WHERE ";
                    for (auto p: conf) {
                        query << spacer << p.key << "='" << p.value.str() << "'";
                        spacer = " AND ";
                    }
                    query << ";";
                    SQLite3Result result = db.query(query.str());

                    if (result.size() == 1) {
                        uuid = result.first()->str();

                        /*query.str(std::string());
                            query << "UPDATE cache_two SET accessed = CURRENT_TIMESTAMP WHERE uuid = '" << uuid << "';";
                            db.exec(query.str());*/ // TODO

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
                        db.exec(query.str());
                    }

                    // check whether .mat and .json file exists and compare settings in program with settings in .json file
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

                    // create .json file if "is_existing" is false
                    if (!is_existing) {
                        conf.save_to_json(path_json.string());
                    }

                    // save everything
                    std::shared_ptr<Hamiltonianmatrix> mat = std::make_shared<Hamiltonianmatrix>(totalmatrix);
                    mat->addFilename(path_mat.string());
                    mat->addIsExisting(is_existing);
                    params.push_back(std::make_shared<Configuration>(conf)); // TODO
                    matrix_path.push_back(path.string());
                    matrix_step.push_back(step_two);
                    matrix_blocks.push_back(blocks_mat_single.size());
                    matrix_block.push_back(idx_block);
                    matrix_dimension.push_back(mat->num_basisvectors());
                    matrix.push_back(std::move(mat));

                    std::cout << "Two-atom Hamiltonian, " <<  matrix.size() << ". Hamiltonian assembled" << std::endl;
                }
            }

            std::cout << ">>TOT" << std::setw(7) << matrix.size() << std::endl; // TODO

            std::cout << "Two-atom Hamiltonian, all Hamiltonians assembled" << std::endl;
        }

        // --- diagonalize matrices using MpiLoadbalancingSimple ---
        run(matrix, matrix_diag);
    }

    std::shared_ptr<const BasisnamesTwo> names() const{
        return basis_two;
    }

private:
    std::shared_ptr<const HamiltonianOne> hamiltonian_one1;
    std::shared_ptr<const HamiltonianOne> hamiltonian_one2;
    std::shared_ptr<BasisnamesTwo> basis_two;
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


// ############################################################################
// ### MAIN LOOP ##############################################################
// ############################################################################

int main(int argc, char **argv) {
    std::cout << std::unitbuf;

    // === Parse command line ===
    boost::filesystem::path path_config;
    boost::filesystem::path path_cache;
    int c;
    opterr = 0;
    while((c = getopt (argc, argv, "c:o:")) != -1) {
        switch (c) {
        case 'c':
            path_config = boost::filesystem::absolute(optarg);
            break;
        case 'o':
            path_cache = boost::filesystem::absolute(optarg);
            break;
        case '?':
            if (optopt == 'c') {
                std::cout << "Option \"-"<< static_cast<char>(optopt) <<"\" requires an argument." << std::endl;
            } else {
                std::cout << "Unknown option \"-"<< static_cast<char>(optopt) <<"\"." << std::endl;
            }
            return 1;
        default:
            abort();
        }
    }

    if (not path_config.string().size()) {
        std::cout << "Required option \"-c filename.json\"." << std::endl;
        return 1;
    }

    if (not boost::filesystem::exists(path_config)) {
        std::cout << "Path " << path_config << " does not exist." << std::endl;
        return 1;
    }

    if (not path_cache.string().size()) {
        std::cout << "Required option \"-o filename.json\"." << std::endl;
        return 1;
    }

    if (not boost::filesystem::exists(path_cache)) {
        std::cout << "Path " << path_cache << " does not exist." << std::endl;
        return 1;
    }

    // === Load configuration ===
    Configuration config;
    config.load_from_json(path_config.string());

    bool existAtom1 = config.count("species1") && config.count("n1") && config.count("l1") && config.count("j1") && config.count("m1");
    bool existAtom2 = config.count("species2") && config.count("n2") && config.count("l2") && config.count("j2") && config.count("m2");

    // === Solve the system ===
    auto mpi = std::make_shared<MpiEnvironment>(argc, argv);

    bool combined = config["samebasis"].str() == "true";

    if (combined) {
        if (config["species1"].str() != config["species2"].str()) {
            std::cout << "species1 and species2 has to be the same in order to use the same basis set." << std::endl;
            return 1;
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one;
        if (existAtom1 && existAtom2) {
            if (mpi->rank() == 0) std::cout << ">>TYP" << std::setw(7) << 3 << std::endl;
            auto basisnames_one = std::make_shared<BasisnamesOne>(BasisnamesOne::fromBoth(config));
            hamiltonian_one = std::make_shared<HamiltonianOne>(config, path_cache, mpi, basisnames_one);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && config.count("minR")) {
            if (mpi->rank() == 0) std::cout << ">>TYP" << std::setw(7) << 2 << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, mpi, hamiltonian_one);
        }
    } else {
        std::shared_ptr<HamiltonianOne> hamiltonian_one1;
        if (existAtom1) {
            if (mpi->rank() == 0) std::cout << ">>TYP" << std::setw(7) << 0 << std::endl;
            auto basisnames_one1 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromFirst(config));
            hamiltonian_one1 = std::make_shared<HamiltonianOne>(config, path_cache, mpi, basisnames_one1);
        }
        std::shared_ptr<HamiltonianOne> hamiltonian_one2;
        if (existAtom2) {
            if (mpi->rank() == 0) std::cout << ">>TYP" << std::setw(7) << 1 << std::endl;
            auto basisnames_one2 = std::make_shared<BasisnamesOne>(BasisnamesOne::fromSecond(config));
            hamiltonian_one2 = std::make_shared<HamiltonianOne>(config, path_cache, mpi, basisnames_one2);
        }
        std::shared_ptr<HamiltonianTwo> hamiltonian_two;
        if (existAtom1 && existAtom2 && config.count("minR")) {
            if (mpi->rank() == 0) std::cout << ">>TYP" << std::setw(7) << 2 << std::endl;
            hamiltonian_two = std::make_shared<HamiltonianTwo>(config, path_cache, mpi, hamiltonian_one1, hamiltonian_one2);
        }
    }

    // hamiltonian_two->saveLines(); //TODO

    // === Communicate that everything is finished ===
    MPI_Barrier(mpi->world());
    if (mpi->rank() == 0) std::cout << ">>END" << std::endl;


    //StateTwo startstate(config);

    //HamiltonianOne hamiltonian_one(mpi, config, startstate);


    //HamiltonianOne hamiltonian_one1(mpi, config, startstate.first(), startstate.second());
    /*HamiltonianTwo hamiltonian_two(mpi, config, hamiltonian_one1, hamiltonian_one2);

    if (mpi->rank() == 0) {
        if (hamiltonian_two.names().size() < 20) {
            std::cout << std::endl ;

            Eigen::IOFormat CleanFmt(2, 0, "  ", "\n", "", "");

            auto mat = hamiltonian_two.get(0);
            std::cout << Eigen::MatrixXd(mat->entries()).format(CleanFmt) << std::endl<< std::endl;
            std::cout << Eigen::MatrixXd(mat->basis()).format(CleanFmt) << std::endl<< std::endl;

            mat = hamiltonian_two.get(1);
            std::cout << Eigen::MatrixXd(mat->entries()).format(CleanFmt) << std::endl<< std::endl;
            std::cout << Eigen::MatrixXd(mat->basis()).format(CleanFmt) << std::endl<< std::endl;

            for (const auto &state : hamiltonian_two.names()) {
                std::cout << state << std::endl;
            }
        }
    }*/

    return 0;
}
