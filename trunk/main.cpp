#include "dtypes.h"
#include "MpiEnvironment.h"
#include "MpiLoadbalancingComplex.h"
#include "MpiLoadbalancingSimple.h"
#include "Vectorizable.h"
#include "Serializable.h"
#include "DipoleMatrix.hpp"
#include "QuantumDefect.hpp"
#include "Basisnames.h"
#include "MatrixElements.h"

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

#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>


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

    Hamiltonianmatrix(size_t szBasis, size_t szCoordinates) : Serializable()  {
        triplets_basis.reserve(szBasis);
        triplets_entries.reserve(szCoordinates);
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

    Hamiltonianmatrix abs() const {
        return Hamiltonianmatrix(entries_.cwiseAbs().cast<scalar_t>(), basis_);
    }
    Hamiltonianmatrix changeBasis(eigen_sparse_t basis) const{
        auto transformator = basis_.transpose()*basis;
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

    void findUnnecessaryStates(std::vector<bool> &isNecessary) const {
        for (eigen_idx_t k=0; k<basis_.outerSize(); ++k) {
            for (eigen_iterator_t it(basis_,k); it; ++it) {
                isNecessary[it.row()] = true;
            }
        }
    }

    void removeUnnecessaryStates(const std::vector<bool> &isNecessary) {
        bytes.clear();

        // build transformator
        std::vector<eigen_triplet_t> triplets_transformator;
        triplets_transformator.reserve(num_coordinates());

        size_t idxCoordinate = 0;
        for (size_t idx = 0; idx < this->num_coordinates(); ++idx) {
            if (isNecessary[idx]) {
                triplets_transformator.push_back(eigen_triplet_t(idxCoordinate++,idx,1));
            }
        }

        eigen_sparse_t transformator(idxCoordinate,this->num_coordinates());
        transformator.setFromTriplets(triplets_transformator.begin(), triplets_transformator.end());

        // apply transformator
        basis_ = transformator*basis_;
    }

    friend Hamiltonianmatrix combineSym(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs) {
        return lhs.duplicate(SYM, rhs);
    }
    friend Hamiltonianmatrix combineAsym(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs) {
        return lhs.duplicate(ASYM, rhs);
    }
    friend Hamiltonianmatrix combineAll(const Hamiltonianmatrix &lhs, const Hamiltonianmatrix &rhs) {
        return lhs.duplicate(ALL, rhs);
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
            std::cout << "The data type used in the program does not fit the data type used in the serialized objects." << std::endl;
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
    void setFilename(std::string name) {
        filename_ = name;
    }

    std::string& filename() {
        return filename_;
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
    void addSymetrized(mode_t mode, std::vector<idx_t> mapping, idx_t row_1, idx_t row_2, idx_t col_1, idx_t col_2, scalar_t val, std::vector<eigen_triplet_t> &triplets_entries) const {
        idx_t row = mapping[this->num_basisvectors()*row_1 + row_2];
        idx_t col = mapping[this->num_basisvectors()*col_1 + col_2];
        if((mode == ALL) || (mode == SYM && row_1 <= row_2 && col_1 <= col_2) || (mode == ASYM && row_1 < row_2 && col_1 < col_2)) {
            real_t factor = 1;
            if (mode == SYM && row_1 == row_2) factor *= 1./sqrt(2.);
            if (mode == SYM && col_1 == col_2) factor *= 1./sqrt(2.);
            triplets_entries.push_back(eigen_triplet_t(row, col, factor*val));
        }
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
                    addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);

                    // <1a 2a|I x V|1b 2b> = <2a 1a|V x I|2b 1b>
                    row_1 = unitmatrix_idx;
                    row_2 = hamiltonian_triple.row();
                    col_1 = unitmatrix_idx;
                    col_2 = hamiltonian_triple.col();
                    addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);

                    // --- mixed terms ---
                    if (mode == ALL) continue;

                    // <1a 2a|V x I|2b 1b> = <2a 1a|I x V|1b 2b>
                    row_1 = hamiltonian_triple.row();
                    row_2 = unitmatrix_idx;
                    col_1 = unitmatrix_idx;
                    col_2 = hamiltonian_triple.col();
                    addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);

                    // <1a 2a|I x V|2b 1b> = <2a 1a|V x I|1b 2b>
                    row_1 = unitmatrix_idx;
                    row_2 = hamiltonian_triple.row();
                    col_1 = hamiltonian_triple.col();
                    col_2 = unitmatrix_idx;
                    addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, triplets_entries);
                }
            }
        }

        mat.entries().setFromTriplets(triplets_entries.begin(), triplets_entries.end());

        return mat;
    }

    eigen_sparse_t entries_;
    eigen_sparse_t basis_;
    enum mode_t {ALL, SYM, ASYM};

    bytes_t bytes;
    size_t idxStart;

    std::string filename_;

    std::vector<eigen_triplet_t> triplets_basis;
    std::vector<eigen_triplet_t> triplets_entries;
};

class Hamiltonian : protected MpiLoadbalancingSimple<Hamiltonianmatrix, Hamiltonianmatrix> {
public:
    Hamiltonian(std::shared_ptr<MpiEnvironment> mpi) : MpiLoadbalancingSimple(mpi, 1000) {}
    std::shared_ptr<Hamiltonianmatrix> get(size_t idx) {
        return matrix_diag[idx];
    }
    std::shared_ptr<const Hamiltonianmatrix> get(size_t idx) const{
        return matrix_diag[idx];
    }
    size_t size() const {
        return matrix_diag.size();
    }

protected:
    std::shared_ptr<Hamiltonianmatrix> doProcessing(std::shared_ptr<Hamiltonianmatrix> work) {

        // if results can not be loaded
        if (true or !work->load()) {
            // diagonalization
            Eigen::SelfAdjointEigenSolver<eigen_dense_t> eigensolver(eigen_dense_t(work->entries()));

            // eigenvalues
            eigen_vector_t evals = eigensolver.eigenvalues().cast<scalar_t>();
            work->entries().setZero();
            work->entries().reserve(evals.size());
            for (eigen_idx_t idx = 0; idx < evals.size(); ++idx) {
                work->entries().insert(idx, idx) = evals.coeffRef(idx);
            }
            work->entries().makeCompressed();

            // eigenvectors
            //work->basis() = work->basis()*eigensolver.eigenvectors().sparseView();
            //work->basis().prune(1e-4,0.5);
            work->basis() = (work->basis()*eigensolver.eigenvectors().sparseView(1e-4,0.5)).pruned(1e-4,0.5); //TODO

            // save result
            work->save();

            std::cout << evals.size() << std::endl;
        }
        return work;
    }

    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix;
    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix_diag;
};

class HamiltonianOne : public Hamiltonian{
public:
    HamiltonianOne(std::shared_ptr<MpiEnvironment> mpi, const StateOne &startstate1 , const StateOne &startstate2) : Hamiltonian(mpi), basis_one(startstate1, startstate2) {
        build();
    }

    HamiltonianOne(std::shared_ptr<MpiEnvironment> mpi, const StateOne &startstate) : Hamiltonian(mpi), basis_one(startstate) {
        build();
    }

    const BasisnamesOne& names() const {
        return basis_one;
    }

protected:

    template <class T>
    void changeToSpherical(real_t val_x, real_t val_y, real_t val_z, T& val_p, T& val_m, T& val_0) {
        std::cout << "test" << val_x << std::endl;
        if(val_x != 0 || val_y != 0) {
            std::cout << "For fields with non-zero x,y-coordinates, a complex data type is needed." << std::endl;
            abort();
        }
        val_p = val_x;
        val_m = val_y;
        val_0 = val_z;
    }

    template <class T>
    void changeToSpherical(real_t val_x, real_t val_y, real_t val_z, std::complex<T>& val_p, std::complex<T>& val_m, std::complex<T>& val_0) {
        val_p = std::complex<T>(-val_x/std::sqrt(2),val_y/std::sqrt(2));
        val_m = std::complex<T>(val_x/std::sqrt(2),val_y/std::sqrt(2));
        val_0 = std::complex<T>(val_z,0);
    }

    void build() {
        if (mpi->rank() == 0) {
            // if not distant dependent, nSteps should be 1 here
            size_t nSteps = 7*4;

            real_t tol = 1e-32;

            real_t energycutoff = 0.7e-5;

            real_t min_E_x = 0;
            real_t min_E_y = 0;
            real_t min_E_z = 0;
            real_t max_E_x = 0;
            real_t max_E_y = 0;
            real_t max_E_z = 1e-11;
            real_t min_B_x = 0;
            real_t min_B_y = 0;
            real_t min_B_z = 50*4.254382e-10;
            real_t max_B_x = 0;
            real_t max_B_y = 0;
            real_t max_B_z = 50*4.254382e-10;

            // === calculate one-atom Hamiltonians ===
            scalar_t min_E_0, min_E_p, min_E_m, min_B_0, min_B_p, min_B_m, max_E_0, max_E_p, max_E_m, max_B_0, max_B_p, max_B_m;
            changeToSpherical(min_E_x, min_E_y, min_E_z, min_E_p, min_E_m, min_E_0);
            changeToSpherical(max_E_x, max_E_y, max_E_z, max_E_p, max_E_m, max_E_0);
            changeToSpherical(min_B_x, min_B_y, min_B_z, min_B_p, min_B_m, min_B_0);
            changeToSpherical(max_B_x, max_B_y, max_B_z, max_B_p, max_B_m, max_B_0);

            bool exist_E_0 = (std::abs(min_E_0) != 0 || std::abs(max_E_0) != 0);
            bool exist_E_p = (std::abs(min_E_p) != 0 || std::abs(max_E_p) != 0);
            bool exist_E_m = (std::abs(min_E_m) != 0 || std::abs(max_E_m) != 0);
            bool exist_B_0 = (std::abs(min_B_0) != 0 || std::abs(max_B_0) != 0);
            bool exist_B_p = (std::abs(min_B_p) != 0 || std::abs(max_B_p) != 0);
            bool exist_B_m = (std::abs(min_B_m) != 0 || std::abs(max_B_m) != 0);

            // --- count entries of Hamiltonian parts ---
            size_t size_basis = basis_one.size();
            size_t size_energy = basis_one.size();
            size_t size_d_0 = 0;
            size_t size_d_p = 0;
            size_t size_d_m = 0;
            size_t size_S_0 = 0;
            size_t size_S_p = 0;
            size_t size_S_m = 0;
            size_t size_L_0 = 0;
            size_t size_L_p = 0;
            size_t size_L_m = 0;

            std::cout << 8.1 << std::endl;

            for (const auto &state_col : basis_one) {
                for (const auto &state_row : basis_one) {
                    if (state_row.idx < state_col.idx) {
                        continue;
                    }

                    if (exist_E_0 && selectionRulesDipole(state_row, state_col, 0) ) {
                        size_d_0++;
                    } else if (exist_E_p && selectionRulesDipole(state_row, state_col, 1) ) {
                        size_d_p++;
                    } else if (exist_E_m && selectionRulesDipole(state_row, state_col, -1) ) {
                        size_d_m++;
                    }

                    if (exist_B_0 && selectionRulesMomentum(state_row, state_col, 0) ) {
                        size_S_0++;
                        size_L_0++;
                    } else if (exist_B_p && selectionRulesMomentum(state_row, state_col, 1) ) {
                        size_S_p++;
                        size_L_p++;
                    } else if (exist_B_m && selectionRulesMomentum(state_row, state_col, -1) ) {
                        size_S_m++;
                        size_L_m++;
                    }
                }
            }

            std::cout << 8.2 << std::endl;

            // --- construct energy Hamiltonian part ---
            Hamiltonianmatrix hamiltonian_energy(size_basis, size_energy);

            real_t energy_initial = 0;
            for (const auto &state: basis_one.initial()) {
                energy_initial += energy_level("Rb",state.n,state.l,state.j);
            }
            energy_initial /= basis_one.initial().size();

            std::vector<bool> is_necessary(basis_one.size(),false);
            idx_t idx = 0;
            for (const auto &state : basis_one) {
                real_t val = energy_level("Rb",state.n,state.l,state.j)-energy_initial;
                if (std::abs(val) <= energycutoff) {
                    is_necessary[state.idx] = true;
                    hamiltonian_energy.addEntries(idx,idx,val);
                    hamiltonian_energy.addBasis(idx,idx,1);
                    ++idx;
                }
            }
            basis_one.removeUnnecessaryStates(is_necessary);

            hamiltonian_energy.compress(basis_one.dim(), basis_one.dim());

            std::cout << basis_one.size() << std::endl;

            // --- precalculate matrix elements ---
            MatrixElements matrix_elements;

            if (exist_E_0 || exist_E_p || exist_E_m || exist_B_0 || exist_B_p || exist_B_m) {
                matrix_elements.precalculate(basis_one, exist_E_0, exist_E_p, exist_E_m, exist_B_0, exist_B_p, exist_B_m);
            }

            // --- construct field Hamiltonian parts ---
            Hamiltonianmatrix hamiltonian_d_0(size_basis, size_d_0);
            Hamiltonianmatrix hamiltonian_d_p(size_basis, size_d_p);
            Hamiltonianmatrix hamiltonian_d_m(size_basis, size_d_m);
            Hamiltonianmatrix hamiltonian_m_0(size_basis, size_S_0);
            Hamiltonianmatrix hamiltonian_m_p(size_basis, size_S_p);
            Hamiltonianmatrix hamiltonian_m_m(size_basis, size_S_m);

            for (const auto &state_col : basis_one) {
                for (const auto &state_row : basis_one) {
                    if (state_row.idx < state_col.idx) {
                        continue;
                    }

                    if (state_row.idx == state_col.idx) {
                        hamiltonian_d_0.addBasis(state_row.idx,state_col.idx,1);
                        hamiltonian_d_p.addBasis(state_row.idx,state_col.idx,1);
                        hamiltonian_d_m.addBasis(state_row.idx,state_col.idx,1);
                        hamiltonian_m_0.addBasis(state_row.idx,state_col.idx,1);
                        hamiltonian_m_p.addBasis(state_row.idx,state_col.idx,1);
                        hamiltonian_m_m.addBasis(state_row.idx,state_col.idx,1);
                    }

                    if (exist_E_0 && selectionRulesDipole(state_row, state_col, 0) ) {
                        real_t val = matrix_elements.getDipole(state_row, state_col);
                        if (std::abs(val) > tol) {
                            hamiltonian_d_0.addEntries(state_row.idx,state_col.idx,val);
                        }
                    } else if (exist_E_p && selectionRulesDipole(state_row, state_col, 1) ) {
                        real_t val = matrix_elements.getDipole(state_row, state_col);
                        if (std::abs(val) > tol) {
                            hamiltonian_d_p.addEntries(state_row.idx,state_col.idx,val);
                        }
                    } else if (exist_E_m && selectionRulesDipole(state_row, state_col, -1) ) {
                        real_t val = matrix_elements.getDipole(state_row, state_col);
                        if (std::abs(val) > tol) {
                            hamiltonian_d_m.addEntries(state_row.idx,state_col.idx,val);
                        }
                    }

                    if (exist_B_0 && selectionRulesMomentum(state_row, state_col, 0) ) {
                        real_t val = matrix_elements.getMomentum(state_row, state_col);
                        if (std::abs(val) > tol) {
                            hamiltonian_m_0.addEntries(state_row.idx,state_col.idx,val);
                        }
                    } else if (exist_B_p && selectionRulesMomentum(state_row, state_col, 1) ) {
                        real_t val = matrix_elements.getMomentum(state_row, state_col);
                        if (std::abs(val) > tol) {
                            hamiltonian_m_p.addEntries(state_row.idx,state_col.idx,val);
                        }
                    } else if (exist_B_m && selectionRulesMomentum(state_row, state_col, -1) ) {
                        real_t val = matrix_elements.getMomentum(state_row, state_col);
                        if (std::abs(val) > tol) {
                            hamiltonian_m_m.addEntries(state_row.idx,state_col.idx,val);
                        }
                    }
                }
            }

            std::cout << 8.4 << std::endl;

            hamiltonian_d_0.compress(basis_one.dim(), basis_one.dim());
            hamiltonian_d_p.compress(basis_one.dim(), basis_one.dim());
            hamiltonian_d_m.compress(basis_one.dim(), basis_one.dim());
            hamiltonian_m_0.compress(basis_one.dim(), basis_one.dim());
            hamiltonian_m_p.compress(basis_one.dim(), basis_one.dim());
            hamiltonian_m_m.compress(basis_one.dim(), basis_one.dim());

            std::cout << 8.5 << std::endl;

            // --- construct Hamiltonians ---
            matrix.reserve(nSteps);

            for (size_t step = 0; step < nSteps; ++step) {
                real_t normalized_position = step/(nSteps-1.);
                auto mat = std::make_shared<Hamiltonianmatrix>(hamiltonian_energy
                                                               -hamiltonian_d_0*(min_E_0+normalized_position*(max_E_0-min_E_0))
                                                               -hamiltonian_d_p*(min_E_p+normalized_position*(max_E_p-min_E_p))
                                                               -hamiltonian_d_m*(min_E_m+normalized_position*(max_E_m-min_E_m))
                                                               +hamiltonian_m_0*(min_B_0+normalized_position*(max_B_0-min_B_0))
                                                               +hamiltonian_m_p*(min_B_p+normalized_position*(max_B_p-min_B_p))
                                                               +hamiltonian_m_m*(min_B_m+normalized_position*(max_B_m-min_B_m))
                                                               );

                std::stringstream s;
                s << "output/hamiltonian_one_" << step << ".mat";
                mat->setFilename(s.str());

                matrix.push_back(std::move(mat));
            }

            std::cout << 8.6 << std::endl;
        }

        // === diagonalize matrices using MpiLoadbalancingSimple ===
        run(matrix, matrix_diag);













        /*std::vector<eigen_triplet_t> triplets_energy;
            std::vector<eigen_triplet_t> triplets_d_0;
            std::vector<eigen_triplet_t> triplets_d_p;
            std::vector<eigen_triplet_t> triplets_d_m;
            std::vector<eigen_triplet_t> triplets_basis;
            triplets_energy.reserve(basis_one.size()*basis_one.size()*0.1);
            triplets_d_0.reserve(basis_one.size()*basis_one.size()*0.1);
            triplets_d_p.reserve(basis_one.size()*basis_one.size()*0.1);
            triplets_d_m.reserve(basis_one.size()*basis_one.size()*0.1);
            triplets_basis.reserve(basis_one.size());

            std::cout << 2.5 << std::endl;

            // loop over basis states
            for (const auto &state_col : basis_one) {
                for (const auto &state_row : basis_one) {
                    if (state_row.idx < state_col.idx) {
                        continue;
                    }

                    // add entries
                    if (state_row.idx == state_col.idx) {
                        real_t val = energy_level("Rb",state_row.n,state_row.l,state_row.j);
                        triplets_energy.push_back(eigen_triplet_t(state_row.idx,state_col.idx,val));
                    } else if (selection_rules(state_row.l, state_row.j, state_row.m, state_col.l, state_col.j, state_col.m) ) {
                        int order = 1;
                        int q_pol = 0;
                        real_t val = radial_element("Rb", state_row.n, state_row.l, state_row.j, order, "Rb", state_col.n, state_col.l, state_col.j) *
                                angular_element(state_row.l, state_row.j, state_row.m, state_col.l, state_col.j, state_col.m, q_pol);

                        if (fabs(val) > 1e-32) { // TODO
                            triplets_d_0.push_back(eigen_triplet_t(state_row.idx,state_col.idx,val));
                        }
                    }

                    // add basis
                    if (state_row.idx == state_col.idx) {
                        triplets_basis.push_back(eigen_triplet_t(state_row.idx,state_col.idx,1));
                    }
                }
            }

            std::cout << 2.6 << std::endl;

            Hamiltonianmatrix hamiltonian_energy(basis_one.dim(),basis_one.dim());
            Hamiltonianmatrix hamiltonian_d_0(basis_one.dim(),basis_one.dim());

            hamiltonian_energy.entries().setFromTriplets(triplets_energy.begin(), triplets_energy.end());
            hamiltonian_energy.basis().setFromTriplets(triplets_basis.begin(), triplets_basis.end());
            hamiltonian_d_0.entries().setFromTriplets(triplets_d_0.begin(), triplets_d_0.end());
            hamiltonian_d_0.basis().setFromTriplets(triplets_basis.begin(), triplets_basis.end());

            std::cout << 2.7 << std::endl;

            // loop over distances
            for (size_t step = 0; step < nSteps; ++step) {
                auto mat = std::make_shared<Hamiltonianmatrix>(hamiltonian_energy+hamiltonian_d_0*step*1e-11);

                std::stringstream s;
                s << "output/hamiltonian_one_" << step << ".mat";
                mat->setFilename(s.str());

                matrix.push_back(std::move(mat));
            }
            std::cout << 2.8 << std::endl;
        }*/

    }

private:
    BasisnamesOne basis_one;
};

class HamiltonianTwo : public Hamiltonian {
public:
    HamiltonianTwo(std::shared_ptr<MpiEnvironment> mpi, const HamiltonianOne &hamiltonian_one1, const HamiltonianOne &hamiltonian_one2) : Hamiltonian(mpi),  basis_two(hamiltonian_one1.names(), hamiltonian_one2.names()) {
        assert(hamiltonian_one1.size() == hamiltonian_one2.size());

        if (mpi->rank() == 0) {
            bool distant_dependent = (hamiltonian_one1.size() == 1) ? false : true;
            size_t nSteps = (distant_dependent) ? hamiltonian_one1.size() : 100; // TODO

            real_t energycutoff = 5; // TODO

            std::cout << 1 << std::endl;

            // --- load one-atom Hamiltonians ---
            std::vector<Hamiltonianmatrix> free_sym_list;
            std::vector<Hamiltonianmatrix> free_asym_list;
            free_sym_list.reserve(hamiltonian_one1.size());
            free_asym_list.reserve(hamiltonian_one1.size());

            std::vector<bool> isNecessary(basis_two.size(),false);

            for (size_t i = 0; i < hamiltonian_one1.size(); ++i) {
                free_sym_list.push_back(combineSym(*hamiltonian_one1.get(i), *hamiltonian_one2.get(i)));
                free_asym_list.push_back(combineAsym(*hamiltonian_one1.get(i), *hamiltonian_one2.get(i)));

                // use energy cutoff
                free_sym_list.back().applyCutoff(energycutoff);
                free_asym_list.back().applyCutoff(energycutoff);

                // find states that became unnecessary
                free_sym_list.back().findUnnecessaryStates(isNecessary);
                free_asym_list.back().findUnnecessaryStates(isNecessary);
            }

            // remove unnecessary states
            basis_two.removeUnnecessaryStates(isNecessary);

            std::cout << 2 << std::endl;

            // --- calculate two-atom Hamiltonians ---
            std::vector<eigen_triplet_t> triplets_entries_dipdip;
            std::vector<eigen_triplet_t> triplets_entries_dipquad;
            std::vector<eigen_triplet_t> triplets_entries_quadquad;
            std::vector<eigen_triplet_t> triplets_basis_dipdip;
            std::vector<eigen_triplet_t> triplets_basis_dipquad;
            std::vector<eigen_triplet_t> triplets_basis_quadquad;
            triplets_entries_dipdip.reserve(basis_two.size()*basis_two.size());
            triplets_entries_dipquad.reserve(basis_two.size()*basis_two.size());
            triplets_entries_quadquad.reserve(basis_two.size()*basis_two.size());
            triplets_basis_dipdip.reserve(basis_two.size());
            triplets_basis_dipquad.reserve(basis_two.size());
            triplets_basis_quadquad.reserve(basis_two.size());

            std::cout << 2.5 << std::endl;

            // loop over basis states
            for (const auto &state_col : basis_two) {
                for (const auto &state_row : basis_two) {

                    // add entries
                    if (true) { // check for selection rules // TODO
                        real_t val = (rand() % 100 - 50)/10.; // calculate value of matrix element // TODO
                        triplets_entries_dipdip.push_back(eigen_triplet_t(state_row.idx,state_col.idx,val));
                    }
                    if (true) { // check for selection rules // TODO
                        real_t val = 0; // calculate value of matrix element // TODO
                        triplets_entries_dipquad.push_back(eigen_triplet_t(state_row.idx,state_col.idx,val));
                    }
                    if (true) { // check for selection rules // TODO
                        real_t val = 0; // calculate value of matrix element // TODO
                        triplets_entries_quadquad.push_back(eigen_triplet_t(state_row.idx,state_col.idx,val));
                    }

                    // add basis
                    if (state_row.idx == state_col.idx) {
                        triplets_basis_dipdip.push_back(eigen_triplet_t(state_row.idx,state_col.idx,1));
                        triplets_basis_dipquad.push_back(eigen_triplet_t(state_row.idx,state_col.idx,1));
                        triplets_basis_quadquad.push_back(eigen_triplet_t(state_row.idx,state_col.idx,1));
                    }
                }
            }

            std::cout << 2.6 << std::endl;

            Hamiltonianmatrix interaction_dipdip(basis_two.dim(),basis_two.dim());
            Hamiltonianmatrix interaction_dipquad(basis_two.dim(),basis_two.dim());
            Hamiltonianmatrix interaction_quadquad(basis_two.dim(),basis_two.dim());

            std::cout << 2.7 << std::endl;

            interaction_dipdip.entries().setFromTriplets(triplets_entries_dipdip.begin(), triplets_entries_dipdip.end());
            interaction_dipdip.basis().setFromTriplets(triplets_basis_dipdip.begin(), triplets_basis_dipdip.end());
            interaction_dipquad.entries().setFromTriplets(triplets_entries_dipquad.begin(), triplets_entries_dipquad.end());
            interaction_dipquad.basis().setFromTriplets(triplets_basis_dipquad.begin(), triplets_basis_dipquad.end());
            interaction_quadquad.entries().setFromTriplets(triplets_entries_quadquad.begin(), triplets_entries_quadquad.end());
            interaction_quadquad.basis().setFromTriplets(triplets_basis_quadquad.begin(), triplets_basis_quadquad.end());

            std::cout << 3 << std::endl;

            auto interaction_dipdip_sym = interaction_dipdip.changeBasis(free_sym_list[0].basis());
            auto interaction_dipquad_sym = interaction_dipquad.changeBasis(free_sym_list[0].basis());
            auto interaction_quadquad_sym = interaction_quadquad.changeBasis(free_sym_list[0].basis());
            auto interaction_dipdip_asym = interaction_dipdip.changeBasis(free_asym_list[0].basis());
            auto interaction_dipquad_asym = interaction_dipquad.changeBasis(free_asym_list[0].basis());
            auto interaction_quadquad_asym = interaction_quadquad.changeBasis(free_asym_list[0].basis());

            std::cout << 4 << std::endl;

            // --- search for submatrices ---
            auto test_sym = interaction_dipdip_sym.abs()+interaction_dipquad_sym.abs()+interaction_quadquad_sym.abs();
            auto test_asym = interaction_dipdip_asym.abs()+interaction_dipquad_asym.abs()+interaction_quadquad_asym.abs();
            for (size_t i = 0; i < hamiltonian_one1.size(); ++i) {
                test_sym += free_sym_list[i].abs();
                test_asym += free_asym_list[i].abs();
            }

            (void) test_sym; // TODO
            (void) test_asym; // TODO

            // neglect irrelevant submatrices
            // TODO

            size_t nSubmatrices = 2; // TODO

            // save submatrices (treat sym/asym as submatrices, too)
            std::vector<Hamiltonianmatrix> arr_interaction_dipdip;
            std::vector<Hamiltonianmatrix> arr_interaction_dipquad;
            std::vector<Hamiltonianmatrix> arr_interaction_quadquad;
            arr_interaction_dipdip.reserve(nSubmatrices);
            arr_interaction_dipquad.reserve(nSubmatrices);
            arr_interaction_quadquad.reserve(nSubmatrices);

            arr_interaction_dipdip.push_back(interaction_dipdip_sym);
            arr_interaction_dipquad.push_back(interaction_dipquad_sym);
            arr_interaction_quadquad.push_back(interaction_quadquad_sym);
            arr_interaction_dipdip.push_back(interaction_dipdip_asym);
            arr_interaction_dipquad.push_back(interaction_dipquad_asym);
            arr_interaction_quadquad.push_back(interaction_quadquad_asym);

            std::vector<std::vector<Hamiltonianmatrix>> arr_free_list(hamiltonian_one1.size());
            for (size_t i = 0; i < hamiltonian_one1.size(); ++i) {
                arr_free_list[i].reserve(nSubmatrices);

                arr_free_list[i].push_back(free_sym_list[i]);
                arr_free_list[i].push_back(free_asym_list[i]);
            }

            std::cout << 5 << std::endl;

            // --- construct total Hamiltonians ---
            matrix.reserve(nSteps*nSubmatrices);

            size_t i = 0;

            // loop over distances
            for (size_t step = 0; step < nSteps; ++step) {
                real_t distance = 1*step/30.+1; // TODO

                // loop over submatrices
                for (size_t sub = 0; sub < nSubmatrices; ++sub) {

                    // add the two-atom Hamiltonians to the one-atom Hamiltonians
                    auto mat = std::make_shared<Hamiltonianmatrix>(arr_free_list[i][sub] +
                                                                   1./pow(distance,3.)*arr_interaction_dipdip[sub] +
                                                                   1./pow(distance,4.)*arr_interaction_dipquad[sub] +
                                                                   1./pow(distance,5.)*arr_interaction_quadquad[sub]);

                    std::stringstream s;
                    s << "output/hamiltonian_two_" << step << "_" << sub << ".mat";
                    mat->setFilename(s.str());

                    // store the total Hamiltonian
                    matrix.push_back(std::move(mat));
                }

                // select new one-atom Hamiltonians, if they are distant dependent
                if (distant_dependent && i+1 < hamiltonian_one1.size()) ++i;
            }

            std::cout << 6 << std::endl;
        }

        // --- diagonalize matrices using MpiLoadbalancingSimple ---
        run(matrix, matrix_diag);
    }

    const BasisnamesTwo& names() const{
        return basis_two;
    }

private:
    BasisnamesTwo basis_two;
};


// ############################################################################
// ### MAIN LOOP ##############################################################
// ############################################################################

int main(int argc, char **argv) {
    std::cout << std::unitbuf;

    auto mpi = std::make_shared<MpiEnvironment>(argc, argv);

    StateTwo startstate({{66,66}}, {{0,0}}, {{0.5,0.5}}, {{0.5,0.5}}, {{0.5,0.5}}); // n, l, s, j, m // TODO

    HamiltonianOne hamiltonian_one2(mpi, startstate.second());
    //HamiltonianOne hamiltonian_one1(mpi, startstate.first(), startstate.second());
    /*HamiltonianTwo hamiltonian_two(mpi, hamiltonian_one1, hamiltonian_one2);

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
