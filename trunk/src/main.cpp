#include "dtypes.h"
#include "MpiEnvironment.h"
#include "MpiLoadbalancingComplex.h"
#include "MpiLoadbalancingSimple.h"
#include "MatrixCOO.h"
#include "MatrixCRS.h"
#include "Vectorizable.h"
#include "Serializable.h"

#include <memory>
#include <tuple>
#include <algorithm>
#include <iterator>

#include <iostream>
#include <vector>
#include <math.h>

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


///////////////////// Stuff to implement  /////////////////////

class Eigensystem : private MpiLoadbalancingSimple {
    VecHamiltonian diaghamiltonians_;

public:
    Eigensystem(std::shared_ptr<MpiEnvironment> mpi, VecHamiltonian &hamiltonians) : MpiLoadbalancingSimple(mpi, 1000), diaghamiltonians_(hamiltonians.num_positions(), hamiltonians.num_blocks()) {
        std::vector<std::shared_ptr<Serializable>> vecIn;
        std::vector<std::shared_ptr<Serializable>> vecOut;

        if (mpi->rank() == 0) {
            vecIn->reserve(hamiltonians.size());
            vecOut->reserve(hamiltonians.size());

            for (auto &p: hamiltonians) {
                vecIn.push_back(p);
            }

            for (auto &p: diaghamiltonians_) {
                vecOut.push_back(p);
            }
        }

        auto bufferSerializable = std::make_shared<Hamiltonian>();

        run(vecIn, vecOut, bufferSerializable);
    }

    std::vector<std::vector<real_t>> potentials() {
        std::vector<std::vector<idx_t>> potentials;
        for (size_t n = 0; n<hamiltonians.num_blocks(); ++n) {
            for (auto &p: diaghamiltonians_.GetBlockbasis(n)) {
                // TODO connect basis elements with largest overlapp
            }
            // TODO save to potentials, how the basis elements are connected
        }
        return potentials;
    }

private:
    std::shared_ptr<Serializable> doProcessing(std::shared_ptr<Serializable> work) {

        auto hamiltonian = std::static_pointer_cast<Hamiltonian>(work);
        auto filepath = hamiltonian->filepath();

        if (hamiltonian->is_buffered()) {
            return std::make_shared<Hamiltonian>(filepath);

        } else {
            auto matrix = hamiltonian->matrix();
            auto basis = hamiltonian->basis();

            MatrixCRS diag;
            MatrixCRS evecs;
            diagonalize(matrix, diag, evecs);

            auto diaghamiltonian = std::make_shared<Hamiltonian>(diag, basis*evecs, filepath);
            diaghamiltonian->save();
            return diaghamiltonian;
        }
    }

    void diagonalize(MatrixCRS &matrix, MatrixCRS &diag, MatrixCRS &evecs) {
        // TODO implement the diagonalization
    }
};
*/


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class DemoSimple : private MpiLoadbalancingSimple<MatrixCRS, MatrixCRS> { // <TIn, TOut>
    std::vector<std::shared_ptr<MatrixCRS>> diaghamiltonians_;

public:
    DemoSimple(std::shared_ptr<MpiEnvironment> mpi, std::vector<std::shared_ptr<MatrixCRS>> &hamiltonians) : MpiLoadbalancingSimple(mpi, 1000) {
        run(hamiltonians, diaghamiltonians_);
    }

    std::vector<std::shared_ptr<MatrixCRS>>& diaghamiltonians() {
        return diaghamiltonians_;
    }

private:
    std::shared_ptr<MatrixCRS> doProcessing(std::shared_ptr<MatrixCRS> work) { // zeroth slave of slave group
        work->multiplyScalar(-1.);
        return work;
    }
};


class Demo : private MpiLoadbalancingComplex<MatrixCRS, MatrixCRS, MatrixCOO> { // <TIn, TOut, TVec>
    std::vector<std::shared_ptr<MatrixCRS>> diaghamiltonians_;

public:
    Demo(std::shared_ptr<MpiEnvironment> mpi, std::vector<std::shared_ptr<MatrixCRS>> &hamiltonians) : MpiLoadbalancingComplex(mpi, 10) {
        run(hamiltonians, diaghamiltonians_);
    }

    std::vector<std::shared_ptr<MatrixCRS>>& diaghamiltonians() {
        return diaghamiltonians_;
    }

private:
    std::shared_ptr<MatrixCOO> doMainprocessing(std::shared_ptr<MatrixCRS> work) { // all slaves within slave group
        work->multiplyScalar(-2.);
        return std::make_shared<MatrixCOO>(work->toCOO());
    }

    std::shared_ptr<MatrixCRS> doPostprocessing(std::shared_ptr<MatrixCOO> resultCombined, std::shared_ptr<MatrixCRS> work) { // zeroth slave of slave group
        resultCombined->setDimensions(work->getDimensions());
        return std::make_shared<MatrixCRS>(resultCombined->toCRS());
    }
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include <array>


class State {
public:
    State(int n_, int l_, float s_, float j_, float m_) {
        n[0] = n_;
        l[0] = l_;
        s[0] = s_;
        j[0] = j_;
        m[0] = m_;
        n[1] = 0;
        l[1] = 0;
        s[1] = 0;
        j[1] = 0;
        m[1] = 0;
    }
    State(std::array<int, 2> n_, std::array<int, 2> l_, std::array<float, 2> s_, std::array<float, 2> j_, std::array<float, 2> m_) {
        n = n_;
        l = l_;
        s = s_;
        j = j_;
        m = m_;
    }
protected:
    std::array<int, 2> n, l;
    std::array<float, 2> s, j, m;
};


// ----------------------------------------

class Basisnames {
public:
    Basisnames() {}
    size_t size() {
        return size_;
    }
    State& get(size_t idx) {
        return names[idx];
    }

protected:
    size_t size_;
    std::vector<State> names;
};


class BasisnamesOne : public Basisnames{
public:
    BasisnamesOne() {
        size_ = 4;
        names.reserve(size_);

        for (size_t i=0; i < size_; ++i) { // loop over quantum numbers
            names.push_back(State(i,i,i,i,i));
        }
    }
};


class BasisnamesTwo : public Basisnames{
public:
    BasisnamesTwo(BasisnamesOne basis_one) {
        size_ = basis_one.size()*basis_one.size();
        names.reserve(size_);

        for (size_t idx_1 = 0; idx_1 < basis_one.size(); ++idx_1) {
            auto state_1 = basis_one.get(idx_1);

            for (size_t idx_2 = 0; idx_2 < basis_one.size(); ++idx_2) {
                auto state_2 = basis_one.get(idx_2);

                (void) state_1;
                (void) state_2;

                names.push_back(State(1,1,1,1,1));
            }
        }
    }
};

// ----------------------------------------

class Hamiltonianmatrix : public Serializable {
public:
    Hamiltonianmatrix() : Serializable() {}
    Hamiltonianmatrix(size_t nBasis, size_t nCoordinates, size_t size) : Serializable(), num_basisvectors_(nBasis), num_coordinates_(nCoordinates), entries_(nBasis,nBasis,size), basis_(nCoordinates,nBasis,size)  {}
    Hamiltonianmatrix(MatrixCRS entries, MatrixCRS basis) : Serializable(), num_basisvectors_(basis.getNumCols()), num_coordinates_(basis.getNumRows()), entries_(entries), basis_(basis) {}
    MatrixCRS& entries() {
        return entries_;
    }
    MatrixCRS& basis() {
        return basis_;
    }
    std::array<MatrixCRS,2> basis_symmetric() {
        size_t num_coordinates_single = sqrt(num_coordinates_);

        // --- construct vector of basis indices ---
        std::vector<std::array<size_t,2>> basisindices;
        basisindices.reserve(num_coordinates_);

        for (size_t idx_1 = 0; idx_1 < num_coordinates_single; ++idx_1) {
            for (size_t idx_2 = 0; idx_2 < num_coordinates_single; ++idx_2) {
                basisindices.push_back({idx_1,idx_2});
            }
        }

        // --- construct symmetrization matrix ---
        MatrixCRS sym_crs((num_coordinates_+num_coordinates_single)/2,num_coordinates_,num_coordinates_);
        MatrixCRS asym_crs((num_coordinates_-num_coordinates_single)/2,num_coordinates_,num_coordinates_-num_coordinates_single);

        idx_t sym_idx_row = 0;
        idx_t asym_idx_row = 0;

        for (auto basisindices_row : basisindices) {

            if (basisindices_row[0] > basisindices_row[1]) {
                continue;
            }

            short parts = 0;

            for (auto basisindices_col : basisindices) {
                idx_t idx_col = num_coordinates_single*basisindices_col[0] + basisindices_col[1];

                if (basisindices_row[0] == basisindices_row[1] && basisindices_row[0] == basisindices_col[0] && basisindices_row[1] == basisindices_col[1]) { // e.g 11, 11
                    sym_crs.add(sym_idx_row, idx_col, 1);

                    sym_idx_row++;
                    break;

                } else if (basisindices_row[0] == basisindices_col[0] && basisindices_row[1] == basisindices_col[1]) { // e.g 12, 12
                    sym_crs.add(sym_idx_row, idx_col, 1./sqrt(2.));
                    asym_crs.add(asym_idx_row, idx_col, 1./sqrt(2.));

                    parts++;

                } else if (basisindices_row[0] == basisindices_col[1] && basisindices_row[1] == basisindices_col[0]) { // e.g 12, 21
                    sym_crs.add(sym_idx_row, idx_col, 1./sqrt(2.));
                    asym_crs.add(asym_idx_row, idx_col, -1./sqrt(2.));

                    parts++;
                }

                if (parts == 2) {
                    sym_idx_row++;
                    asym_idx_row++;
                    break;
                }
            }
        }

        return std::array<MatrixCRS,2>({sym_crs,asym_crs});
    }
    void add(Hamiltonianmatrix matrix) {
        (void) matrix; // assimilate basis and add another matrix
    }
    void addSymetrized(mode_t mode, std::vector<idx_t> mapping, idx_t row_1, idx_t row_2, idx_t col_1, idx_t col_2, real_t val, MatrixCOO &coo) {
        idx_t row = mapping[num_basisvectors_*row_1 + row_2];
        idx_t col = mapping[num_basisvectors_*col_1 + col_2];
        if((mode == ALL) || (mode == SYM && row_1 <= row_2 && col_1 <= col_2) || (mode == ASYM && row_1 < row_2 && col_1 < col_2)) {
            double factor = 1;
            if (mode == SYM && row_1 == row_2) factor *= 1./sqrt(2.);
            if (mode == SYM && col_1 == col_2) factor *= 1./sqrt(2.);
            coo.add(row, col, factor*val);
        }
    }
    friend Hamiltonianmatrix operator+(Hamiltonianmatrix lhs, const Hamiltonianmatrix& rhs) {
        return rhs; // TODO
    }
    friend Hamiltonianmatrix operator*(real_t lhs, const Hamiltonianmatrix& rhs) {
        return rhs; // TODO
    }
    Hamiltonianmatrix duplicateSym() {
        return duplicate(SYM);
    }
    Hamiltonianmatrix duplicateAsym() {
        return duplicate(ASYM);
    }
    Hamiltonianmatrix duplicateAll() {
        return duplicate(ALL);
    }
    Hamiltonianmatrix duplicate(mode_t mode) {
        // --- mapping ---
        idx_t i = 0;
        std::vector<idx_t> mapping(num_basisvectors_*num_basisvectors_);
        for (idx_t idx_1 = 0; idx_1 < num_basisvectors_; ++idx_1) {
            for (idx_t idx_2 = 0; idx_2 < num_basisvectors_; ++idx_2) {
                if ((mode != ALL && idx_1 < idx_2) || (mode == ALL) || (mode == SYM && idx_1 == idx_2)) {
                    idx_t idx = num_coordinates_*idx_1 + idx_2;
                    mapping[idx] = i++;
                }
            }
        }

        // --- duplicate basis_ ---
        MatrixCOO basis_coo(num_coordinates_*num_coordinates_,i,basis_.size()*basis_.size()); //TODO

        for (auto triple_1 : basis_) {
            for (auto triple_2 : basis_) {
                if ((mode != ALL && triple_1.col < triple_2.col)) {
                    idx_t idx_row1 = num_basisvectors_*triple_1.row + triple_2.row; // coord1
                    idx_t idx_row2 = num_basisvectors_*triple_2.row + triple_1.row; // coord2
                    idx_t idx_col = mapping[num_coordinates_*triple_1.col + triple_2.col]; // vec

                    int factor = (mode == ASYM) ? -1 : 1;
                    basis_coo.add(idx_row1, idx_col, triple_1.val*triple_2.val/sqrt(2.));
                    basis_coo.add(idx_row2, idx_col, triple_1.val*triple_2.val/sqrt(2.)*factor);

                } else if ((mode == ALL) || (mode == SYM && triple_1.col == triple_2.col)) {
                    idx_t idx_row = num_basisvectors_*triple_1.row + triple_2.row; // coord
                    idx_t idx_col = mapping[num_coordinates_*triple_1.col + triple_2.col]; // vec
                    basis_coo.add(idx_row, idx_col, triple_1.val*triple_2.val);
                }
            }
        }

        basis_coo.sumup();

        // --- duplicate entries_ ---
        MatrixCOO entries_coo(i,i,2*entries_.size()*num_basisvectors_); //TODO

        for (auto hamiltonian_triple : entries_) {
            for (size_t unitmatrix_idx = 0; unitmatrix_idx < num_basisvectors_; ++unitmatrix_idx) {
                idx_t row_1, row_2, col_1, col_2;
                real_t val = hamiltonian_triple.val;

                // --- ordered terms ---
                // <1a 2a|V x I|1b 2b> = <2a 1a|I x V|2b 1b>
                row_1 = hamiltonian_triple.row;
                row_2 = unitmatrix_idx;
                col_1 = hamiltonian_triple.col;
                col_2 = unitmatrix_idx;
                addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, entries_coo);

                // <1a 2a|I x V|1b 2b> = <2a 1a|V x I|2b 1b>
                row_1 = unitmatrix_idx;
                row_2 = hamiltonian_triple.row;
                col_1 = unitmatrix_idx;
                col_2 = hamiltonian_triple.col;
                addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, entries_coo);

                // --- mixed terms ---
                if (mode == ALL) continue;

                // <1a 2a|V x I|2b 1b> = <2a 1a|I x V|1b 2b>
                row_1 = hamiltonian_triple.row;
                row_2 = unitmatrix_idx;
                col_1 = unitmatrix_idx;
                col_2 = hamiltonian_triple.col;
                addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, entries_coo);

                // <1a 2a|I x V|2b 1b> = <2a 1a|V x I|1b 2b>
                row_1 = unitmatrix_idx;
                row_2 = hamiltonian_triple.row;
                col_1 = hamiltonian_triple.col;
                col_2 = unitmatrix_idx;
                addSymetrized(mode, mapping, row_1, row_2, col_1, col_2, val, entries_coo);
            }
        }

        entries_coo.sumup();

        return Hamiltonianmatrix(entries_coo.toCRS(), basis_coo.toCRS());
    }
    Hamiltonianmatrix changeBasis(MatrixCRS basis) {
        auto entries = entries_.transformBackward(basis_).transformForward(basis);
        return Hamiltonianmatrix(entries, basis);
    }

    bytes_t serialize() {
        bytes_t bytes_entries = entries_.serialize();
        bytes_t bytes_basis = basis_.serialize();
        size_t entries_size = bytes_entries.size();
        size_t basis_size = bytes_basis.size();

        bytes_t bytes(2*sizeof(size_t));

        auto pbytes = bytes.begin();
        serializeItem(pbytes,entries_size);
        serializeItem(pbytes,basis_size);

        bytes.reserve(bytes.size() + bytes_entries.size() + bytes_basis.size() );
        bytes.insert(bytes.end(), bytes_entries.begin(), bytes_entries.end());
        bytes.insert(bytes.end(), bytes_basis.begin(), bytes_basis.end());

        return bytes;
    }
    void deserialize(bytes_t &bytes) {
        size_t entries_size;
        size_t basis_size;

        auto pbytes = bytes.begin();
        deserializeItem(pbytes,entries_size);
        deserializeItem(pbytes,basis_size);

        bytes_t bytes_entries(pbytes, pbytes+entries_size);
        bytes_t bytes_basis(pbytes+entries_size, pbytes+entries_size+basis_size);

        entries_.deserialize(bytes_entries);
        basis_.deserialize(bytes_basis);

        num_basisvectors_ = basis_.getNumCols();
        num_coordinates_ = basis_.getNumRows();
    }

protected:
    size_t num_basisvectors_;
    size_t num_coordinates_;
    MatrixCRS entries_;
    MatrixCRS basis_;
    enum mode_t {ALL, SYM, ASYM};


};

class Hamiltonian : protected MpiLoadbalancingSimple<Hamiltonianmatrix, Hamiltonianmatrix> {
public:
    Hamiltonian(std::shared_ptr<MpiEnvironment> mpi) : MpiLoadbalancingSimple(mpi, 1000) {}
    std::shared_ptr<Hamiltonianmatrix> get(size_t idx) {
        return matrix_diag[idx];
    }
    size_t size() const {
        return matrix_diag.size();
    }

protected:
    std::shared_ptr<Hamiltonianmatrix> doProcessing(std::shared_ptr<Hamiltonianmatrix> work) {
        work->entries().multiplyScalar(-1.); // Diagonalization
        return work;
    }

    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix;
    std::vector<std::shared_ptr<Hamiltonianmatrix>> matrix_diag;
};

class HamiltonianOne : public Hamiltonian{
public:
    HamiltonianOne(std::shared_ptr<MpiEnvironment> mpi, BasisnamesOne basis_one) : Hamiltonian(mpi) {

        if (mpi->rank() == 0) {
            // if not distant dependent, nSteps should be 1 here
            size_t nSteps = 10;
            matrix.reserve(nSteps);

            // loop over distances
            for (size_t idx_i = 0; idx_i < nSteps; ++idx_i) {
                auto mat = std::make_shared<Hamiltonianmatrix>(basis_one.size(),basis_one.size(),basis_one.size()*basis_one.size()*0.5);

                for (size_t idx_row = 0; idx_row < basis_one.size(); ++idx_row) {
                    auto state_row = basis_one.get(idx_row);

                    for (size_t idx_col = 0; idx_col < basis_one.size(); ++idx_col) {
                        auto state_col = basis_one.get(idx_col);

                        (void) state_row;
                        (void) state_col;

                        // add entries
                        if (idx_row == idx_col) { // check for selection rules // TODO
                            real_t val = idx_row; // calculate value of matrix element // TODO
                            mat->entries().add(idx_row,idx_col,val);
                        }

                        // add basis
                        if (idx_row == idx_col) {
                            mat->basis().add(idx_row,idx_col,1);
                        }
                    }
                }

                matrix.push_back(std::move(mat));
            }
        }

        // --- diagonalize matrices using MpiLoadbalancingSimple ---
        run(matrix, matrix_diag);
    }
};

class HamiltonianTwo : public Hamiltonian {
public:
    HamiltonianTwo(std::shared_ptr<MpiEnvironment> mpi, BasisnamesTwo basis_two, HamiltonianOne hamiltonian_one) : Hamiltonian(mpi) {

        if (mpi->rank() == 0) {
            bool distant_dependent = (hamiltonian_one.size() == 1) ? false : true;
            size_t nSteps = (distant_dependent) ? hamiltonian_one.size() : 10;

            // --- load one-atom Hamiltonians ---
            auto free_sym = hamiltonian_one.get(0)->duplicateSym();
            auto free_asym = hamiltonian_one.get(0)->duplicateAsym();

            // use energy cuttoff
            // TODO

            // --- calculate two-atom Hamiltonians ---
            Hamiltonianmatrix interaction_dipdip(basis_two.size(),basis_two.size(),25);
            Hamiltonianmatrix interaction_dipquad(basis_two.size(),basis_two.size(),25);
            Hamiltonianmatrix interaction_quadquad(basis_two.size(),basis_two.size(),25);

            for (idx_t idx_row = 0; idx_row < basis_two.size(); ++idx_row) {
                auto state_row = basis_two.get(idx_row);

                for (idx_t idx_col = 0; idx_col < basis_two.size(); ++idx_col) {
                    auto state_col = basis_two.get(idx_col);

                    (void) state_row;
                    (void) state_col;

                    // add entries
                    if (true) { // check for selection rules // TODO
                        real_t val = 2; // calculate value of matrix element // TODO
                        interaction_dipdip.entries().add(idx_row,idx_col,val);
                    }
                    if (true) { // check for selection rules // TODO
                        real_t val = 2; // calculate value of matrix element // TODO
                        interaction_dipquad.entries().add(idx_row,idx_col,val);
                    }
                    if (true) { // check for selection rules // TODO
                        real_t val = 2; // calculate value of matrix element // TODO
                        interaction_quadquad.entries().add(idx_row,idx_col,val);
                    }

                    // add basis
                    if (idx_row == idx_col) {
                        interaction_dipdip.basis().add(idx_row,idx_col,1);
                        interaction_dipquad.basis().add(idx_row,idx_col,1);
                        interaction_quadquad.basis().add(idx_row,idx_col,1);
                    }
                }
            }

            auto interaction_dipdip_sym = interaction_dipdip.changeBasis(free_sym.basis());
            auto interaction_dipquad_sym = interaction_dipquad.changeBasis(free_sym.basis());
            auto interaction_quadquad_sym = interaction_quadquad.changeBasis(free_sym.basis());

            auto interaction_dipdip_asym = interaction_dipdip.changeBasis(free_asym.basis());
            auto interaction_dipquad_asym = interaction_dipquad.changeBasis(free_asym.basis());
            auto interaction_quadquad_asym = interaction_quadquad.changeBasis(free_asym.basis());

            // --- search for submatrices ---
            auto test_sym = free_sym+interaction_dipdip_sym+interaction_dipquad_sym+interaction_quadquad_sym; // TODO Watch out for distant dependent free_sym/free_asym
            auto test_asym = free_asym+interaction_dipdip_asym+interaction_dipquad_asym+interaction_quadquad_asym;
            (void) test_sym; // TODO
            (void) test_asym; // TODO

            // neglect irrelevant submatrices
            // TODO

            size_t nSubmatrices = 2; // TODO

            // save submatrices (treat sym/asym as submatrices, too)
            std::vector<Hamiltonianmatrix> arr_free;
            std::vector<Hamiltonianmatrix> arr_interaction_dipdip;
            std::vector<Hamiltonianmatrix> arr_interaction_dipquad;
            std::vector<Hamiltonianmatrix> arr_interaction_quadquad;
            arr_free.reserve(nSubmatrices);
            arr_interaction_dipdip.reserve(nSubmatrices);
            arr_interaction_dipquad.reserve(nSubmatrices);
            arr_interaction_quadquad.reserve(nSubmatrices);

            arr_free.push_back(free_sym); // TODO
            arr_interaction_dipdip.push_back(interaction_dipdip_sym);
            arr_interaction_dipquad.push_back(interaction_dipquad_sym);
            arr_interaction_quadquad.push_back(interaction_quadquad_sym);

            arr_free.push_back(free_asym); // TODO
            arr_interaction_dipdip.push_back(interaction_dipdip_asym);
            arr_interaction_dipquad.push_back(interaction_dipquad_asym);
            arr_interaction_quadquad.push_back(interaction_quadquad_asym);

            // --- construct total Hamiltonians ---
            matrix.reserve(nSteps*nSubmatrices);

            // loop over distances
            for (size_t step = 0; step < nSteps; ++step) {
                real_t distance = 1; // TODO

                // loop over submatrices
                for (size_t sub = 0; sub < nSubmatrices; ++sub) {

                    // add the two-atom Hamiltonians to the one-atom Hamiltonians
                    auto mat = arr_free[sub] + 1./pow(distance,3.)*arr_interaction_dipdip[sub] + 1./pow(distance,4.)*arr_interaction_dipquad[sub] + 1./pow(distance,5.)*arr_interaction_quadquad[sub];

                    // store the total Hamiltonian
                    matrix.push_back(std::move(std::make_shared<Hamiltonianmatrix>(mat)));
                }

                // load new one-atom Hamiltonians, if they are distant dependent
                if (distant_dependent && step+1 < hamiltonian_one.size()) {
                    free_sym = hamiltonian_one.get(step+1)->duplicateSym();
                    free_asym = hamiltonian_one.get(step+1)->duplicateAsym();
                }
            }
        }

        // --- diagonalize matrices using MpiLoadbalancingSimple ---
        run(matrix, matrix_diag);

        // --- output ---
        if (mpi->rank() == 0) {
            for (auto m: matrix_diag) {
                m->basis().print();
                std::cout << "------------" << std::endl;
            }
        }
    }
};


// ----------------------------------------










// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// ############################################################################
// ### MAIN LOOP ##############################################################
// ############################################################################

int main(int argc, char **argv) {
    auto mpi = std::make_shared<MpiEnvironment>(argc, argv);

    BasisnamesOne basis_one;
    HamiltonianOne hamiltonian_one(mpi, basis_one);

    BasisnamesTwo basis_two(basis_one);
    HamiltonianTwo hamiltonian_two(mpi, basis_two, hamiltonian_one);


    //BasisnamesOne basis_two(basis_one);


    /*std::vector<std::shared_ptr<MatrixCRS>> hamiltonians;
    if (mpi->rank() == 0) {
        for (int i = 0; i < 3; ++i) {
            auto crs = std::make_shared<MatrixCRS>(5,4,5);
            crs->add(1,0,2);
            crs->add(1,0,3);
            crs->add(0,2,4);
            crs->add(0,2,2);
            crs->add(0,1,10*i);
            crs->sumup();
            hamiltonians.push_back(std::move(crs));
        }
    }

    auto demo_simple = DemoSimple(mpi, hamiltonians);
    auto demo = Demo(mpi, hamiltonians);

    if (mpi->rank() == 0) {
        std::cout << std::endl << std::endl;
        for(auto &p: demo_simple.diaghamiltonians()) {
            p->print();
            std::cout <<  "--------------------" << std::endl;
        }
        std::cout << std::endl << std::endl;

        for(auto &p: demo.diaghamiltonians()) {
            p->print();
            std::cout <<  "--------------------" << std::endl;
        }
        std::cout << std::endl << std::endl;
    }


    if (mpi->rank() == 0) {
        MatrixCOO coo(4,4,5);
        coo.add(2,1,1);
        coo.add(1,0,2);
        coo.add(2,0,3);
        coo.add(2,2,4);
        coo.add(2,2,2);
        coo.sumup();
        coo.print();
        std::cout << "-------" << std::endl;

        auto crs = coo.toCRS();
        crs->print();
        std::cout << "-------" << std::endl;

        auto coo2 = crs->toCOO();
        coo2->print();
        std::cout << "-------" << std::endl;

        MatrixCRS crs2(4,4,5);
        crs2.add(1,0,2);
        crs2.add(1,0,3);
        crs2.add(0,2,4);
        crs2.add(0,2,2);
        crs2.add(0,1,1);
        crs2.sumup();
        crs2.print();
        std::cout << "-------" << std::endl;

        auto tmp = crs2.serialize();
        MatrixCRS crs3(tmp);
        crs3.print();
        std::cout << "-------" << std::endl;

    }*/

    return 0;
}
