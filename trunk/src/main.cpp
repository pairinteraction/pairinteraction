#define _MPI_CPP_BINDINGS
#include <mpi.h>
#include <cassert>
#include <memory>
#include <tuple>
#include <algorithm>
#include <iterator>

#include <iostream>
#include <vector>
#include <math.h>


#define REAL MPI_FLOAT
typedef float real;
typedef unsigned int idx;

enum {WORKTAG = 0, DIETAG = 1};

namespace definitions {
int numSlaves = 2;
}

// ############################################################################
// ### FUNCTIONS ##############################################################
// ############################################################################

// ////////////////////////////////////////////////////////////////////////////
// /// Master /////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////

int masterIO(MPI_Comm comm, std::vector<int> vecSlave0) {
    std::vector<std::vector<real>> vecWork(100,std::vector<real>(5));
    std::vector<real> collection;
    int msgSize;
    MPI_Status status;
    MPI_Request request;
    unsigned int numWorkingSlaves0 = 0;

    // === Loop over work ===
    for (auto &work: vecWork) {
        if (numWorkingSlaves0 < vecSlave0.size()) {
            // --- Send new work ---
            MPI_Isend(&work[0], work.size(), REAL, vecSlave0[numWorkingSlaves0], WORKTAG, comm, &request);

            // --- Increase number of working slaves0 ---
            numWorkingSlaves0 += 1;
        } else {
            // --- Receive message from slaves ---
            // Probe for an incoming message from master
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

            // Get message size
            MPI_Get_count(&status, MPI_INT, &msgSize);

            // Allocate buffer
            collection.resize(msgSize);

            // Receive message
            MPI_Recv(&collection[0], collection.size(), REAL, status.MPI_SOURCE, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

            // --- Send new work ---
            MPI_Isend(&work[0], work.size(), REAL, status.MPI_SOURCE, WORKTAG, comm, &request);

            // --- Process message ---
            // TODO
        }
    }

    // === Collect - until now unreceived - results ===
    for (unsigned int i = 0; i < numWorkingSlaves0; i++) {
        // --- Receive message from slaves ---
        // Probe for an incoming message from master
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);

        // Get message size
        MPI_Get_count(&status, MPI_INT, &msgSize);

        // Allocate buffer
        collection.resize(msgSize);

        // Receive message
        MPI_Recv(&collection[0], collection.size(), REAL, status.MPI_SOURCE, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

        // --- Process message ---
        // TODO
    }

    // === Kill all slaves0 ===
    for (auto &dest: vecSlave0) {
        MPI_Isend(NULL, 0, MPI_INT, dest, DIETAG, comm, &request);
    }

    return 0;
}

// ////////////////////////////////////////////////////////////////////////////
// /// Slave //////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////

int slaveIO(MPI_Comm comm, int myRank, int mySlave, MPI_Comm myComm){
    MPI_Status status;
    std::vector<real> work, result, collection;
    std::vector<int> displacement;
    int msgSize, mySize;

    MPI_Comm_size(myComm, &mySize);

    // === Do work ===
    while (true) {
        // --- Receive message from master ---
        if (myRank == 0) {
            // Probe for an incoming message from master
            MPI_Probe(0, MPI_ANY_TAG, comm, &status);

            if (status.MPI_TAG != DIETAG) {
                // Get message size
                MPI_Get_count(&status, MPI_INT, &msgSize);

                // Allocate buffer
                work.resize(msgSize);

                // Receive message
                MPI_Recv(&work[0], work.size(), REAL, 0, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
            } else {
                // Set message size to -1 to account for the DIETAG
                msgSize = -1;
            }
        }

        // --- Receive message from slave0 ---
        if (mySize > 1) {
            // Broadcast message size
            MPI_Bcast(&msgSize, 1, MPI_INT, 0, myComm);

            if (msgSize != -1) {
                if ( myRank != 0 ) {
                    // Allocate buffer
                    work.resize(msgSize);
                }

                // Broadcast message
                MPI_Bcast(&work[0], work.size(), REAL, 0, myComm);
            }
        }

        // --- Stop working in case of DIETAG ---
        if (msgSize == -1) {
            break;
        }

        // --- Do calculations ---
        int  rank;
        MPI_Comm_rank(comm, &rank);
        std::cout << "Hello from rank " << rank << ", slave group " << mySlave << ", rank within slave group " << myRank << std::endl;
        result.resize(5);

        // --- Send results to slave0 ---
        if (mySize > 1) {
            // Tell slave0 how many elements one has
            std::vector<int> numElements(mySize);
            int sz = result.size();
            MPI_Gather(&sz, 1, MPI_INT, &numElements[0], 1, MPI_INT, 0, myComm);

            if ( myRank == 0 ) {
                // Calculate displacements
                int nTotal = 0;
                displacement.resize(mySize);
                for (int i = 0; i < mySize; i++) {
                    displacement[i] = nTotal;
                    nTotal += numElements[i];

                }

                // Allocate buffer
                collection.resize(nTotal);
            }

            // Collect everything into slave0
            MPI_Gatherv(&result[0], result.size(), REAL, &collection[0], &numElements[0], &displacement[0], REAL, 0, myComm);
        } else {
            collection = result;
        }

        // --- Send results to master ---
        if (myRank == 0) {
            // Convert results from COO format into CSR format
            // TODO

            // Send results to master
            MPI_Send(&collection[0], collection.size(), REAL, 0, 0, comm);
        }
    }

    return 0;
}


class MpiEnvironment {
    int rank_;
    int size_;
    MPI::Intracomm world_;

public:
    MpiEnvironment(int argc, char **argv) {
        MPI::Init(argc, argv);
        assert(MPI::Is_initialized() == true);
        rank_ = MPI::COMM_WORLD.Get_rank();
        size_ = MPI::COMM_WORLD.Get_size();
        world_ = MPI::COMM_WORLD;
    }

    ~MpiEnvironment() {
        MPI::Finalize();
        assert(MPI::Is_finalized() == true);
    }

    int rank() const {return rank_;}
    int size() const {return size_;}
    const MPI::Intracomm& world() const {return world_;}
};

class Loadbalancing {
    std::shared_ptr<MpiEnvironment> mpi;
    int numSlaves;
    std::vector<int> vecSlave0;

public:
    Loadbalancing(std::shared_ptr<MpiEnvironment> mpi, int num) : mpi(mpi) {
        // === Correct the number of slaves if it is to heigh ===
        numSlaves = std::min(num, mpi->size()-1);

        // === Construct a vector that contains the rank of the 0th slave inside each slave group ===
        int numUnusedProcessors = mpi->size()-1;
        int numUnusedSlaves = numSlaves;
        int idxProcessor = 1;

        myColor = numSlaves;
        vecSlave0.resize(numSlaves);

        // Loop over slaves
        for (int idxSlave = 0; idxSlave<numSlaves; idxSlave++) {
            int numProcessorsSlave = ceil(numUnusedProcessors/numUnusedSlaves);
            numUnusedProcessors -= numProcessorsSlave;
            numUnusedSlaves -= 1;

            vecSlave0[idxSlave] = idxProcessor;

            // Assign same "color" to the slaves that belong together
            if (mpi->rank() >= idxProcessor && mpi->rank()-1 < idxProcessor+numProcessorsSlave) {
                myColor = idxSlave;
            }

            // Increase the index of used processors
            idxProcessor += numProcessorsSlave;
        }

        // === Build slave specific communicator ===
        myWorld = mpi->world().Split(myColor, mpi->rank());
        myRank = myWorld.Get_rank();
        mySize = myWorld.Get_size();

        // === Wait for all slave communicators to be build ===
        mpi->world().Barrier();
    }

private:
    void runMaster(std::vector<std::vector<char>> &dataIn, std::vector<std::vector<char>> &dataOut) {
        std::vector<MPI::Request> requests;
        unsigned int numWorkingSlaves0 = 0;

        dataOut.clear();
        dataOut.reserve(dataIn.size());
        requests.reserve(dataIn.size()*2+std::min(dataIn.size(), vecSlave0.size()));


        // === Loop over work ===
        for (auto &work: dataIn) {
            if (numWorkingSlaves0 < vecSlave0.size()) {
                // --- Send new work ---
                requests.push_back(mpi->world().Isend(&work[0], work.size(), MPI::CHAR, vecSlave0[numWorkingSlaves0], WORKTAG));

                // --- Increase number of working slaves0 ---
                numWorkingSlaves0 += 1;
            } else {
                // --- Receive message from slaves ---

                // Probe for an incoming message from master
                MPI::Status status;
                mpi->world().Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, status);


                // Allocate buffer
                int msgSize = status.Get_count(MPI::CHAR);
                dataOut.push_back(std::vector<char>(msgSize));

                // Receive message
                requests.push_back(mpi->world().Irecv(&dataOut.back(), msgSize, MPI::CHAR, status.Get_source(), MPI::ANY_TAG));

                // --- Send new work ---
                requests.push_back(mpi->world().Isend(&work[0], work.size(), MPI::CHAR, status.Get_source(), WORKTAG));
            }
        }



        // === Collect - until now unreceived - results ===
        for (unsigned int i = 0; i < numWorkingSlaves0; i++) {
            // --- Receive message from slaves ---

            // Probe for an incoming message from master
            MPI::Status status;
            mpi->world().Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, status);

            // Allocate buffer
            int msgSize = status.Get_count(MPI::CHAR);
            dataOut.push_back(std::vector<char>(msgSize));

            // Receive message
            requests.push_back(mpi->world().Irecv(&dataOut.back(), msgSize, MPI::CHAR, status.Get_source(), MPI::ANY_TAG));
        }

        // === Kill all slaves0 ===
        for (auto &dest: vecSlave0) {
            requests.push_back(mpi->world().Isend(NULL, 0, MPI::INT, dest, DIETAG));
        }

        MPI::Request::Waitall(requests.size(), &requests[0]);
    }

    void runSlave() {
        // === Do work ===
        while (true) {
            std::vector<char> work;
            int msgSize;

            // --- Receive message from master ---
            if (myRank == 0) {
                // Probe for an incoming message from master
                MPI::Status status;
                mpi->world().Probe(0, MPI::ANY_TAG, status);

                if (status.Get_tag() != DIETAG) {
                    // Get message size
                    msgSize = status.Get_count(MPI::CHAR);

                    // Allocate buffer
                    work.resize(msgSize);

                    // Receive message
                    mpi->world().Recv(&work[0], work.size(), MPI::CHAR, 0, MPI::ANY_TAG);
                } else {
                    // Receive message
                    mpi->world().Recv(NULL, 0, MPI::CHAR, 0, MPI::ANY_TAG);

                    // Set message size to -1 to account for the DIETAG
                    msgSize = -1;
                }
            }

            // --- Receive message from slave0 ---
            if (mySize > 1) {
                // Broadcast message size
                myWorld.Bcast(&msgSize, 1, MPI::INT, 0);

                if (msgSize != -1) {
                    if ( myRank != 0 ) {
                        // Allocate buffer
                        work.resize(msgSize);
                    }

                    // Broadcast message
                    myWorld.Bcast(&work[0], work.size(), MPI::CHAR, 0);
                }
            }


            // --- Stop working in case of DIETAG ---
            if (msgSize == -1) {

                break;
            }





            // --- Do calculations ---
            std::cout << "Hello from rank " << mpi->rank() << ", slave group " << myColor << ", rank within slave group " << myRank << std::endl;
            auto result = doMainprocessing(work);



            // --- Send results to slave0 ---
            std::vector<char> collection;

            if (mySize > 1) {
                // Tell slave0 how many elements one has
                std::vector<int> numElements(mySize);
                int sz = result.size();
                myWorld.Gather(&sz, 1, MPI::INT, &numElements[0], 1, MPI::INT, 0);

                std::vector<int> displacement(mySize);
                if ( myRank == 0 ) {
                    // Calculate displacements
                    int nTotal = 0;
                    for (int i = 0; i < mySize; i++) {
                        displacement[i] = nTotal;
                        nTotal += numElements[i];
                    }

                    // Allocate buffer
                    collection.resize(nTotal);
                }

                // Collect everything into slave0
                myWorld.Gatherv(&result[0], result.size(), MPI::CHAR, &collection[0], &numElements[0], &displacement[0], MPI::CHAR, 0);
            } else {
                collection = result;
            }

            // --- Send results to master ---
            if (myRank == 0) {
                // Do postprocessing
                auto results = doPostprocessing(collection);

                // Send results to master
                mpi->world().Send(&results[0], results.size(),  MPI::CHAR, 0, 0);
            }


        }
    }

protected:
    int myColor;
    int myRank;
    int mySize;
    MPI::Intracomm myWorld;

    virtual std::vector<char> doMainprocessing(std::vector<char> &work) = 0; // all slaves within slave group
    virtual std::vector<char> doPostprocessing(std::vector<char> &collection) = 0; // zeroth slave of slave group

    void run(std::vector<std::vector<char>> &dataIn, std::vector<std::vector<char>> &dataOut) {
        if (mpi->rank() == 0) {
            this->runMaster(dataIn, dataOut);
        } else {
            this->runSlave();
        }
    }
};

// --------------------------------------------------------------------------------

class MatrixCOO;

class MatrixCRS {
    std::vector<idx> ptr;
    std::vector<idx> col;
    std::vector<real> val;
    idx nRows;
    idx nCols;
    bool ordered;
    bool sumuped;
public:
    MatrixCRS(idx nRows, idx nCols, size_t size);
    MatrixCRS(std::vector<char> bytes);
    void add(idx rIncrement, idx c, real v);
    void order();
    void sumup();
    void print();
    MatrixCOO toCOO();
    std::vector<char> serialize();
};

class MatrixCOO {
    std::vector<std::tuple<idx,idx,real>> triple; // row, col, val
    idx nRows;
    idx nCols;
    bool ordered;
    bool sumuped;
public:
    MatrixCOO(idx nRows, idx nCols, size_t size);
    MatrixCOO(std::vector<char> bytes);
    void add(idx row, idx col, real val);
    void order();
    void sumup();
    void print();
    MatrixCRS toCRS();
    std::vector<char> serialize();
};

MatrixCRS::MatrixCRS(idx nRows, idx nCols, size_t size) :  nRows(nRows), nCols(nCols), ordered(true), sumuped(true) { // empty
    ptr.reserve(nRows+1);
    col.reserve(size);
    val.reserve(size);
    ptr.push_back(0);
}
MatrixCRS::MatrixCRS(std::vector<char> bytes) { // from bytes
    (void) bytes;
}
void MatrixCRS::add(idx rIncrement, idx c, real v) {
    for (idx i = 0; i < rIncrement; ++i) {
        ptr.push_back(col.size());
    }
    col.push_back(c);
    val.push_back(v);
    ordered = false;
    sumuped = false;
}
void MatrixCRS::order() {
    if (!ordered) {
        while(ptr.size() < nRows+1) {
            ptr.push_back(col.size());
        }

        // Sort everything
        std::vector<size_t> indices;
        indices.reserve(col.size());
        for (size_t i = 0; i < col.size(); ++i) {
            indices.push_back(i);
        }

        for (idx r = 0; r < nRows; ++r) {
            idx ptrStart = (r < ptr.size()) ? ptr[r] : col.size();
            idx ptrEnd = (r+1 < ptr.size()) ? ptr[r+1] : col.size();
            std::sort(indices.begin()+ptrStart, indices.begin()+ptrEnd, [this](size_t i1, size_t i2) {return col[i1] < col[i2];});
        }

        auto colTmp(col);
        auto valTmp(val);
        for (size_t i = 0; i < indices.size(); ++i) {
            col[i] = colTmp[indices[i]];
            val[i] = valTmp[indices[i]];

        }

        ordered = true;
    }
}
void MatrixCRS::sumup() {
    this->order();

    if (!sumuped) {

        // Copy
        auto ptrTmp(ptr);
        auto colTmp(col);
        auto valTmp(val);
        ptr.clear();
        col.clear();
        val.clear();
        ptr.reserve(ptrTmp.size());
        col.reserve(colTmp.size());
        val.reserve(valTmp.size());
        ptr.push_back(0);

        // Sum doubled entries
        for (idx r = 0; r < nRows; ++r) {
            idx oldCol = 0;
            real sumVal = 0;

            for (idx p = ((r < ptrTmp.size()) ? ptrTmp[r] : colTmp.size()); p < ((r+1 < ptrTmp.size()) ? ptrTmp[r+1] : colTmp.size()); ++p) {
                idx c = colTmp[p];
                real v = valTmp[p];
                if (oldCol == c) {
                    sumVal += v;
                } else {
                    if (std::abs(sumVal) > 1e-12) {
                        col.push_back(oldCol);
                        val.push_back(sumVal);
                    }
                    sumVal = v;
                    oldCol = c;
                }
            }

            if (std::abs(sumVal) > 1e-12) {
                col.push_back(oldCol);
                val.push_back(sumVal);
            }

            ptr.push_back(col.size());
        }

        sumuped = true;
    }
}
void MatrixCRS::print() {
    this->order();
    this->sumup();

    for (idx r = 0; r < nRows; ++r) {
        idx p = (r < ptr.size()) ? ptr[r] : col.size();
        idx cc = col[p];

        for (idx c = 0; c < nCols; ++c) {
            real v = 0;
            if (p < ((r+1 < ptr.size()) ? ptr[r+1] : col.size()) && cc == c) {
                v = val[p];
                cc = col[++p];
            }
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
}
MatrixCOO MatrixCRS::toCOO() {
    MatrixCOO coo(nRows, nCols, col.size());
    for (idx r = 0; r < nRows; ++r) {
        for (idx p = ((r < ptr.size()) ? ptr[r] : col.size()); p < ((r+1 < ptr.size()) ? ptr[r+1] : col.size()); ++p) {
            idx c = col[p];
            real v = val[p];
            coo.add(r,c,v);
        }
    }
    return coo;
}
std::vector<char> MatrixCRS::serialize() {
    // TODO
    std::vector<char> bytes;
    return bytes;
}

MatrixCOO::MatrixCOO(idx nRows, idx nCols, size_t size) :  nRows(nRows), nCols(nCols), ordered(true), sumuped(true) { // empty
    triple.reserve(size);
}
MatrixCOO::MatrixCOO(std::vector<char> bytes) { // from bytes
    (void) bytes;
}
void MatrixCOO::add(idx row, idx col, real val) {
    triple.push_back(std::make_tuple(row,col,val));
    ordered = false;
    sumuped = false;
}
void MatrixCOO::order() {
    if (!ordered) {
        // Sort everything
        std::sort(triple.begin(), triple.end(),[](std::tuple<idx,idx,real> t1, std::tuple<idx,idx,real> t2) {return std::get<1>(t1) < std::get<1>(t2);});
        std::stable_sort(triple.begin(), triple.end(),[](std::tuple<idx,idx,real> t1, std::tuple<idx,idx,real> t2) {return std::get<0>(t1) < std::get<0>(t2);});

        ordered = true;
    }
}
void MatrixCOO::sumup() {
    this->order();

    if (!sumuped) {

        // Copy
        auto tripleTmp(triple);
        triple.clear();
        triple.reserve(tripleTmp.size());

        // Sum doubled entries
        idx oldRow = 0;
        idx oldCol = 0;
        real sumVal = 0;
        for (auto &t: tripleTmp) {
            idx row;
            idx col;
            real val;
            std::tie(row, col, val) = t;

            if (row == oldRow && col == oldCol) {
                sumVal += val;
            } else {
                if (std::abs(sumVal) > 1e-12) {
                    triple.push_back(std::make_tuple(oldRow,oldCol,sumVal));
                }
                sumVal = val;
                oldRow = row;
                oldCol = col;
            }
        }
        if (std::abs(sumVal) > 1e-12) {
            triple.push_back(std::make_tuple(oldRow,oldCol,sumVal));
        }

        sumuped = true;
    }
}
void MatrixCOO::print() {
    this->order();
    this->sumup();

    idx n = 0;

    idx row;
    idx col;
    real val;
    std::tie(row, col, val) = triple[n];

    for (idx r = 0; r < nRows; ++r) {
        for (idx c = 0; c < nCols; ++c) {
            real v = 0;
            if (row == r && col == c) {
                v = val;
                std::tie(row, col, val) = triple[++n];
            }
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
}
MatrixCRS MatrixCOO::toCRS() {
    std::stable_sort(triple.begin(), triple.end(),[](std::tuple<idx,idx,real> t1, std::tuple<idx,idx,real> t2) {return std::get<0>(t1) < std::get<0>(t2);});

    MatrixCRS crs(nRows, nCols, triple.size());
    idx lastRow = 0;
    for (auto &tuple: triple) {
        idx row;
        idx col;
        real val;
        std::tie (row, col, val) = tuple;
        crs.add(row-lastRow, col, val);
        lastRow = row;
    }
    return crs;
}
std::vector<char> MatrixCOO::serialize() {
    // TODO
    std::vector<char> bytes;
    return bytes;
}

// --------------------------------------------------------------------------------

class Hamiltonian : private Loadbalancing {
    std::shared_ptr<MpiEnvironment> mpi;
    std::vector<std::vector<char>> vecResult;

public:
    Hamiltonian(std::shared_ptr<MpiEnvironment> mpi) : Loadbalancing(mpi, 5) {
        std::vector<std::vector<char>> vecWork(3,std::vector<char>(5));
        this->run(vecWork, vecResult);
    }

private:
    std::vector<char> doMainprocessing(std::vector<char> &work) { // all slaves within slave group
        // e.g. diagonalization, do nothing if buffered (signalized by using a special header in the work data)
        auto resultSingle = work;
        return resultSingle;
    }

    std::vector<char> doPostprocessing(std::vector<char> &resultCombined) { // zeroth slave of slave group
        // e.g. convert results from COO format into CSR format and save results to buffer, load results if buffered (signalized by using a special header in the work data)
        auto result = resultCombined;
        return result;
    }
};


class Eigensystem : private Loadbalancing {
    std::shared_ptr<MpiEnvironment> mpi;
    std::vector<std::vector<char>> vecResult;

public:
    Eigensystem(std::shared_ptr<MpiEnvironment> mpi) : Loadbalancing(mpi, 5) {
        std::vector<std::vector<char>> vecWork(3,std::vector<char>(5));
        this->run(vecWork, vecResult);
    }

private:
    std::vector<char> doMainprocessing(std::vector<char> &work) { // all slaves within slave group
        // e.g. diagonalization, do nothing if buffered (signalized by using a special header in the work data)
        auto resultSingle = work;
        return resultSingle;
    }

    std::vector<char> doPostprocessing(std::vector<char> &resultCombined) { // zeroth slave of slave group
        // e.g. convert results from COO format into CSR format and save results to buffer, load results if buffered (signalized by using a special header in the work data)
        auto result = resultCombined;
        return result;
    }
};

// ############################################################################
// ### MAIN LOOP ##############################################################
// ############################################################################

int main(int argc, char **argv) {

    auto mpi = std::make_shared<MpiEnvironment>(argc, argv);

    auto hamiltonian = std::make_shared<Hamiltonian>(mpi);
    auto eigensystem = std::make_shared<Eigensystem>(mpi);

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

        MatrixCRS crs = coo.toCRS();
        crs.print();
        std::cout << "-------" << std::endl;

        MatrixCOO coo2 = crs.toCOO();
        coo2.print();
        std::cout << "-------" << std::endl;

        MatrixCRS crs2(4,4,5);
        crs2.add(1,0,2);
        crs2.add(1,0,3);
        crs2.add(0,2,4);
        crs2.add(0,2,2);
        crs2.add(0,1,1);
        crs2.sumup();
        crs2.print();
    }

    //std::vector<std::tuple<idx,idx,real>> data = {std::make_tuple(0,2,2),std::make_tuple(2,1,2),std::make_tuple(0,1,2),std::make_tuple(3,1,2),std::make_tuple(7,1,2),std::make_tuple(5,1,2)};

    //hamiltonian->fromCOO(data,3);






    //auto eigensystems = std::make_shared<Eigensystems>(createEigensystems(mpi, hamiltonian)); // TODO add vecR, vecB, vecE

    // TODO http://cpplove.blogspot.co.at/2013/05/my-take-on-serialization-part-ii.html
    // TODO http://www.codeproject.com/Articles/567242/AplusC-b-bplusObjectplusFactory
    // TODO https://stackoverflow.com/questions/5120768/how-to-implement-the-factory-method-pattern-in-c-correctly

    return 0;




    int rank, myRank, mySlave, numProcessorsTotal, numProcessorsSlaves;
    MPI_Comm myComm;
    std::vector<int> vecSlave0;

    // === Initialize MPI and get ranke ===
    MPI_Init(&argc, &argv);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcessorsTotal);
    numProcessorsSlaves = numProcessorsTotal - 1;

    // === Correct the number of slaves if it is to heigh ===
    definitions::numSlaves = std::min(definitions::numSlaves, numProcessorsSlaves);

    // === Construct a vector that contains the rank of the 0th slave inside each slave group ===
    int numUnusedProcessors = numProcessorsSlaves;
    int numUnusedSlaves = definitions::numSlaves;
    int idxProcessor = 1;

    mySlave = definitions::numSlaves;
    vecSlave0.resize(definitions::numSlaves);

    // Loop over slaves
    for (int idxSlave = 0; idxSlave<definitions::numSlaves; idxSlave++) {
        int numProcessorsSlave = ceil(numUnusedProcessors/numUnusedSlaves);
        numUnusedProcessors -= numProcessorsSlave;
        numUnusedSlaves -= 1;

        vecSlave0[idxSlave] = idxProcessor;

        // Assign same "color" to the slaves that belong together
        if (rank >= idxProcessor && rank-1 < idxProcessor+numProcessorsSlave) {
            mySlave = idxSlave;
        }

        // Increase the index of used processors
        idxProcessor += numProcessorsSlave;
    }

    // === Build slave specific communicator ===
    MPI_Comm_split(MPI_COMM_WORLD, mySlave, rank, &myComm);
    MPI_Comm_rank(myComm, &myRank);


    // === Wait for all slave communicators to be build ===
    MPI_Barrier(MPI_COMM_WORLD);



    /*
    // TODOs
    // mpi namespace
    // #define MPI_FLAG    #ifdef MPI_FLAG ... #endif


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

    // Create Eiegnsystem
    Eigensystems eigensystems(hamiltonians, vecR, vecB, vecE); // use loadbalancing WITH slave groups, read and write to cache if necessary
    eigensystems.save(json);

    // Analyze
    potentials = eigensystems.calculatePotentials();
    potentials.save(json);

    overlap = eigensystems.calculateOverlap(superposition);
    overlap.save(json);


    (*) maybe WITH slave groups, too
*/




    // === Start load balancing ===
    if (rank == 0) {
        masterIO(MPI_COMM_WORLD, vecSlave0);
    } else {
        slaveIO(MPI_COMM_WORLD, myRank, mySlave, myComm);
    }

    // === Finalize MPI ===
    MPI_Finalize();

    return 0;
}
