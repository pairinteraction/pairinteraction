#include "dtypes.h"
#include "MpiEnvironment.h"
#include "MpiLoadbalancingComplex.h"
#include "MpiLoadbalancingSimple.h"
#include "Vectorizable.h"
#include "Serializable.h"
#include "DipoleMatrix.hpp"

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


#include "Iter.h"

template<class T> class Basisnames {
public:
    Basisnames() {}

    size_t size() const {
        return names_.size();
    }
    size_t dim() const {
        return names_.back().idx+1;
    }
    T& get (size_t idx) {
        return names_[idx];
    }
    const T& get (size_t idx) const {
        return names_[idx];
    }
    void set (size_t i, const T &v) {
        names_[i] = v;
    }
    Iter<Basisnames, T> begin() const {
        return Iter<Basisnames, T>( this, 0 );
    }
    Iter<Basisnames, T> end() const {
        return Iter<Basisnames, T>( this, names_.size() );
    }

protected:
    std::vector<T> names_;
};


class BasisnamesOne : public Basisnames<StateOne>{
public:
    BasisnamesOne(const StateOne &startstate) {
        size_t size = 4; // TODO
        names_.reserve(size);

        idx_t idx = 0;

        int delta_n = 1; // TODO
        int delta_l = 2; // TODO
        int delta_j = 100; // TODO
        int delta_m = 3; // TODO

        // loop over quantum numbers
        for (int n = fmax(0, startstate.n - delta_n); n <= startstate.n + delta_n; ++n) {
            for (int l = fmax(0, startstate.l - delta_l); l <= fmin(n-1,startstate.l + delta_l); ++l) {
                for (float j = fmax(abs(l - startstate.s), startstate.j - delta_j); j <= fmin(l + startstate.s, startstate.j + delta_j); ++j) {
                    for (float m = fmax(-j, startstate.m - delta_m); m <= fmin(j, startstate.m + delta_m); ++m) {
                        names_.push_back(StateOne(idx++,n,l,startstate.s,j,0));
                    }
                }
            }
        }
        std::cout << names_.size() << std::endl;
        std::cout << "----------------" << std::endl;
    }

    BasisnamesOne(const StateOne &startstate1, const StateOne &startstate2) {
        size_t size = 4; // TODO
        std::unordered_set<StateOne> names_set;
        names_set.reserve(size);

        idx_t idx = 0;

        int delta_n = 1; // TODO
        int delta_l = 2; // TODO
        int delta_j = 100; // TODO
        int delta_m = 3; // TODO

        // loop over quantum numbers of startstate1
        for (int n = fmax(0, startstate1.n - delta_n); n <= startstate1.n + delta_n; ++n) {
            for (int l = fmax(0, startstate1.l - delta_l); l <= fmin(n-1,startstate1.l + delta_l); ++l) {
                for (float j = fmax(abs(l - startstate1.s), startstate1.j - delta_j); j <= fmin(l + startstate1.s, startstate1.j + delta_j); ++j) {
                    for (float m = fmax(-j, startstate1.m - delta_m); m <= fmin(j, startstate1.m + delta_m); ++m) {
                        names_set.insert(StateOne(idx++,n,l,startstate1.s,j,0));
                    }
                }
            }
        }

        // loop over quantum numbers of startstate2
        for (int n = fmax(0, startstate2.n - delta_n); n <= startstate2.n + delta_n; ++n) {
            for (int l = fmax(0, startstate2.l - delta_l); l <= fmin(n-1,startstate2.l + delta_l); ++l) {
                for (float j = fmax(abs(l - startstate2.s), startstate2.j - delta_j); j <= fmin(l + startstate2.s, startstate2.j + delta_j); ++j) {
                    for (float m = fmax(-j, startstate2.m - delta_m); m <= fmin(j, startstate2.m + delta_m); ++m) {
                        names_set.insert(StateOne(idx++,n,l,startstate2.s,j,0));
                    }
                }
            }
        }

        names_ = std::vector<StateOne>(names_set.begin(), names_set.end());
        std::cout << names_.size() << std::endl;
        std::cout << "----------------" << std::endl;
    }
};


class BasisnamesTwo : public Basisnames<StateTwo>{
public:
    BasisnamesTwo(const BasisnamesOne &basis_one1, const BasisnamesOne &basis_one2) {
        size_t size = basis_one1.size()*basis_one2.size();
        names_.reserve(size);

        idx_t idx = 0;

        // loop over single atom states
        for (const auto &state_1 : basis_one1) {
            for (const auto &state_2 : basis_one2) {
                names_.push_back(StateTwo(idx++,state_1,state_2));
            }
        }
    }
    void removeUnnecessaryStates(const std::vector<bool> &isNecessary) {
        auto tmp = names_;
        names_.clear();
        names_.reserve(tmp.size());

        // loop over all two-atom states
        for (auto state : tmp) {
            if (isNecessary[state.idx]) {
                names_.push_back(state);
            }
        }

        names_.shrink_to_fit();
    }
};

// ----------------------------------------

class Hamiltonianmatrix : public Serializable {
public:
    Hamiltonianmatrix() : Serializable() {}
    Hamiltonianmatrix(size_t nBasis, size_t nCoordinates) : Serializable(), entries_(nBasis,nBasis), basis_(nCoordinates,nBasis)  {}
    Hamiltonianmatrix(eigen_sparse_t entries, eigen_sparse_t basis) : Serializable(), entries_(entries), basis_(basis) {}

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

    Hamiltonianmatrix abs() const {
        return Hamiltonianmatrix(entries_.cwiseAbs(), basis_);
    }
    Hamiltonianmatrix changeBasis(eigen_sparse_t basis) const{
        auto transformator = basis_.transpose()*basis;
        auto entries = transformator.transpose()*entries_*transformator;
        return Hamiltonianmatrix(entries, basis);
    }

    void applyCutoff(real_t cutoff) {
        bytes.clear();

        // build transformator
        auto diag = entries_.diagonal();

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
    friend Hamiltonianmatrix operator*(const real_t& lhs,  Hamiltonianmatrix rhs) {
        rhs.bytes.clear();
        rhs.entries_ *= lhs;
        return rhs;
    }
    friend Hamiltonianmatrix operator*(Hamiltonianmatrix lhs,  const real_t& rhs) {
        lhs.bytes.clear();
        lhs.entries_ *= rhs;
        return lhs;
    }
    Hamiltonianmatrix& operator+=(const Hamiltonianmatrix& rhs) {
        bytes.clear();
        entries_ += rhs.entries_;
        return *this;
    }

    bytes_t& serialize() {
        doSerialization();
        return bytes;
    }

    void doSerialization() {
        if (bytes.size() == 0) {
            entries_.makeCompressed();
            basis_.makeCompressed();

            std::vector<char> name(filename_.begin(), filename_.end());
            storage_idx_t entries_mode = (entries_.IsRowMajor) ? 1 : 0; // 0: csc, 1: csr
            storage_idx_t entries_rows = entries_.rows();
            storage_idx_t entries_cols = entries_.cols();
            std::vector<storage_real_t> entries_data(entries_.valuePtr(), entries_.valuePtr()+entries_.nonZeros());
            std::vector<storage_idx_t> entries_indices(entries_.innerIndexPtr(), entries_.innerIndexPtr()+entries_.nonZeros());
            std::vector<storage_idx_t> entries_indptr(entries_.outerIndexPtr(), entries_.outerIndexPtr()+entries_.outerSize());
            storage_idx_t basis_mode = (basis_.IsRowMajor) ? 1 : 0; // 0: csc, 1: csr
            storage_idx_t basis_rows = basis_.rows();
            storage_idx_t basis_cols = basis_.cols();
            std::vector<storage_real_t> basis_data(basis_.valuePtr(), basis_.valuePtr()+basis_.nonZeros());
            std::vector<storage_idx_t> basis_indices(basis_.innerIndexPtr(), basis_.innerIndexPtr()+basis_.nonZeros());
            std::vector<storage_idx_t> basis_indptr(basis_.outerIndexPtr(), basis_.outerIndexPtr()+basis_.outerSize());

            Serializer s;
            s << name;
            idxStart = s.position();
            s << entries_mode;
            s << entries_rows;
            s << entries_cols;
            s << entries_data;
            s << entries_indices;
            s << entries_indptr;
            s << basis_mode;
            s << basis_rows;
            s << basis_cols;
            s << basis_data;
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
        std::vector<char> name;
        storage_idx_t entries_mode, entries_rows, entries_cols;
        std::vector<real_t> entries_data;
        std::vector<idx_t> entries_indices;
        std::vector<idx_t> entries_indptr;
        storage_idx_t basis_mode, basis_rows, basis_cols;
        std::vector<real_t> basis_data;
        std::vector<idx_t> basis_indices;
        std::vector<idx_t> basis_indptr;


        Serializer s;
        s.load(bytes);
        s >> name;
        idxStart = s.position();
        s >> entries_mode;
        s >> entries_rows;
        s >> entries_cols;
        s >> entries_data;
        s >> entries_indices;
        s >> entries_indptr;
        s >> basis_mode;
        s >> basis_rows;
        s >> basis_cols;
        s >> basis_data;
        s >> basis_indices;
        s >> basis_indptr;

        filename_ = std::string(&name[0], name.size());
        entries_ = eigen_sparse_t(entries_rows,entries_cols);
        entries_.makeCompressed();
        entries_.resizeNonZeros(entries_data.size());
        std::copy(entries_data.begin(),entries_data.end(),entries_.valuePtr());
        std::copy(entries_indices.begin(),entries_indices.end(),entries_.innerIndexPtr());
        std::copy(entries_indptr.begin(),entries_indptr.end(),entries_.outerIndexPtr());
        entries_.finalize();
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
    void addSymetrized(mode_t mode, std::vector<idx_t> mapping, idx_t row_1, idx_t row_2, idx_t col_1, idx_t col_2, real_t val, std::vector<eigen_triplet_t> &triplets_entries) const {
        idx_t row = mapping[this->num_basisvectors()*row_1 + row_2];
        idx_t col = mapping[this->num_basisvectors()*col_1 + col_2];
        if((mode == ALL) || (mode == SYM && row_1 <= row_2 && col_1 <= col_2) || (mode == ASYM && row_1 < row_2 && col_1 < col_2)) {
            double factor = 1;
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
                            triplets_basis.push_back(eigen_triplet_t(idx_row1, idx_col, triple_1.value()*triple_2.value()/sqrt(2.)));
                            triplets_basis.push_back(eigen_triplet_t(idx_row2, idx_col, triple_1.value()*triple_2.value()/sqrt(2.)*factor));

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
                    real_t val = hamiltonian_triple.value();

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

        std::cout << work->entries().rows() << std::endl;

        // if results can not be loaded
        if (!work->load()) {

            // diagonalization
            Eigen::SelfAdjointEigenSolver<eigen_dense_t> eigensolver(eigen_dense_t(work->entries()));

            // eigenvalues
            auto evals = eigensolver.eigenvalues();
            work->entries().setZero();
            work->entries().reserve(evals.size());
            for (eigen_idx_t idx = 0; idx < evals.size(); ++idx) {
                work->entries().insert(idx, idx) = evals.coeffRef(idx);
            }
            work->entries().makeCompressed();

            // eigenvectors
            auto tmp = work->basis();
            work->basis() = tmp*eigensolver.eigenvectors();
            work->basis().prune(1e-4,0.5);

            // save result
            work->save();
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
    void build() {
        if (mpi->rank() == 0) {
            // if not distant dependent, nSteps should be 1 here
            size_t nSteps = 1;
            matrix.reserve(nSteps);

            // loop over distances
            for (size_t step = 0; step < nSteps; ++step) {
                std::vector<eigen_triplet_t> triplets_entries;
                std::vector<eigen_triplet_t> triplets_basis;
                triplets_entries.reserve(basis_one.size()*basis_one.size());
                triplets_basis.reserve(basis_one.size());

                // loop over basis states
                for (const auto &state_col : basis_one) {
                    for (const auto &state_row : basis_one) {

                        // add entries
                        if (state_row.idx == state_col.idx) { // check for selection rules // TODO
                            //real_t val = (rand() % 100 - 50)/20.+state_row.idx; // calculate value of matrix element // TODO
                            real_t val = 0;
                            int order = 1;
                            int q_pol = 1;

                            int n1 = state_row.n;
                            int l1 = state_row.l;
                            int j1 = state_row.j;
                            int m1 = state_row.m;

                            int n2 = state_col.n;
                            int l2 = state_col.l;
                            int j2 = state_col.j;
                            int m2 = state_col.m;

                            if ( selection_rules(l1, j1, m1, l2, j2, m2) )
                                val = radial_element("Rb", n1, l1, j1, order, "Rb", n2, l2, j2) * angular_element(l1, j1, m1, l2, j2, m2, q_pol);
                            triplets_entries.push_back(eigen_triplet_t(state_row.idx,state_col.idx,val));
                        }

                        // add basis
                        if (state_row.idx == state_col.idx) {
                            triplets_basis.push_back(eigen_triplet_t(state_row.idx,state_col.idx,1));
                        }
                    }
                }

                auto mat = std::make_shared<Hamiltonianmatrix>(basis_one.dim(),basis_one.dim());
                mat->entries().setFromTriplets(triplets_entries.begin(), triplets_entries.end());
                mat->basis().setFromTriplets(triplets_basis.begin(), triplets_basis.end());

                std::stringstream s;
                s << "output/hamiltonian_one_" << step << ".mat";
                mat->setFilename(s.str());

                matrix.push_back(std::move(mat));
            }
        }

        // --- diagonalize matrices using MpiLoadbalancingSimple ---
        run(matrix, matrix_diag);
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

    StateTwo startstate({120,120}, {4,4}, {0.5,0.5}, {4.5,4.5}, {0.5,0.5}); // n, l, s, j, m // TODO

    HamiltonianOne hamiltonian_one2(mpi, startstate.second());
    HamiltonianOne hamiltonian_one1(mpi, startstate.first(), startstate.second());
    HamiltonianTwo hamiltonian_two(mpi, hamiltonian_one1, hamiltonian_one2);

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
    }

    return 0;
}
