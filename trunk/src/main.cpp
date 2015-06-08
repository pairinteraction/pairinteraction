#include "dtypes.h"
#include "MpiEnvironment.h"
#include "MpiLoadbalancing.h"
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

* mpi namespace
* #define MPI_FLAG    #ifdef MPI_FLAG ... #endif
* use std::vector< std::vector<unsigned char>> instead of std::vector<triple>
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

    std::vector<std::vector<dreal>> potentials() {
        std::vector<std::vector<didx>> potentials;
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

            auto diaghamiltonian = std::make_shared<Hamiltonian>(diag, evecs*basis, filepath); // TODO check, whether calculating evecs*basis is correct
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

class DemoSimple : private MpiLoadbalancingSimple {
    std::vector<std::shared_ptr<MatrixCRS>> diaghamiltonians_;

public:
    DemoSimple(std::shared_ptr<MpiEnvironment> mpi, std::vector<std::shared_ptr<MatrixCRS>> &hamiltonians) : MpiLoadbalancingSimple(mpi, 1000), diaghamiltonians_(hamiltonians) {
        for (auto &p: diaghamiltonians_) p = std::move(std::make_shared<MatrixCRS>());

        std::vector<std::shared_ptr<Serializable>> vecIn (hamiltonians.begin(), hamiltonians.end());
        std::vector<std::shared_ptr<Serializable>> vecOut (diaghamiltonians_.begin(), diaghamiltonians_.end());
        auto bufferSerializable = std::make_shared<MatrixCRS>();

        run(vecIn, vecOut, bufferSerializable);
    }

    std::vector<std::shared_ptr<MatrixCRS>>& diaghamiltonians() {
        return diaghamiltonians_;
    }

private:
    std::shared_ptr<Serializable> doProcessing(std::shared_ptr<Serializable> work) { // zeroth slave of slave group
        auto crs = std::static_pointer_cast<MatrixCRS>(work);
        crs->multiplyScalar(-1.);
        return crs;
    }
};


class Demo : private MpiLoadbalancing {
    std::vector<std::shared_ptr<MatrixCRS>> diaghamiltonians_;

public:
    Demo(std::shared_ptr<MpiEnvironment> mpi, std::vector<std::shared_ptr<MatrixCRS>> &hamiltonians) : MpiLoadbalancing(mpi, 10), diaghamiltonians_(hamiltonians) {
        for (auto &p: diaghamiltonians_) p = std::move(std::make_shared<MatrixCRS>());

        std::vector<std::shared_ptr<Serializable>> vecIn (hamiltonians.begin(), hamiltonians.end());
        std::vector<std::shared_ptr<Serializable>> vecOut (diaghamiltonians_.begin(), diaghamiltonians_.end());
        auto bufferSerializable = std::make_shared<MatrixCRS>();
        auto bufferVectorizable = std::make_shared<MatrixCOO>();

        run(vecIn, vecOut, bufferSerializable, bufferVectorizable);
    }

    std::vector<std::shared_ptr<MatrixCRS>>& diaghamiltonians() {
        return diaghamiltonians_;
    }

private:
    std::shared_ptr<Vectorizable> doMainprocessing(std::shared_ptr<Serializable> work) { // all slaves within slave group
        auto crs = std::static_pointer_cast<MatrixCRS>(work);
        crs->multiplyScalar(-2.);
        return crs->toCOO();
    }

    std::shared_ptr<Serializable> doPostprocessing(std::shared_ptr<Vectorizable> resultCombined, std::shared_ptr<Serializable> work) { // zeroth slave of slave group
        auto coo = std::static_pointer_cast<MatrixCOO>(resultCombined);
        auto crs = std::static_pointer_cast<MatrixCRS>(work);
        coo->setDimensions(crs->getDimensions());
        return coo->toCRS();
    }
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// ############################################################################
// ### MAIN LOOP ##############################################################
// ############################################################################

int main(int argc, char **argv) {
    auto mpi = std::make_shared<MpiEnvironment>(argc, argv);


    std::vector<std::shared_ptr<MatrixCRS>> hamiltonians;
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

    }

    return 0;
}
