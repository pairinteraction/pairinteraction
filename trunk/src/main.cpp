#include <mpi.h>
#include <iostream>
#include <vector>
#include <math.h>


#define REAL MPI_FLOAT
typedef float real;

enum {WORKTAG = 0, DIETAG = 1};

namespace definitions {
   int numSlaves = 2;
}

using namespace std;

// ############################################################################
// ### FUNCTIONS ##############################################################
// ############################################################################

// ////////////////////////////////////////////////////////////////////////////
// /// Master /////////////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////////////////////////

int masterIO(MPI_Comm comm, vector<int> vecSlave0) {
    vector<vector<real>> vecWork(100,vector<real>(5));
    vector<real> collection;
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
    vector<real> work, result, collection;
    vector<int> displacement;
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
        cout << "Hello from rank " << rank << ", slave group " << mySlave << ", rank within slave group " << myRank << endl;
        result.resize(5);

        // --- Send results to slave0 ---
        if (mySize > 1) {
            // Tell slave0 how many elements one has
            vector<int> numElements(mySize);
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

// ############################################################################
// ### MAIN LOOP ##############################################################
// ############################################################################

int main(int argc,char **argv) {
    int rank, myRank, mySlave, numProcessorsTotal, numProcessorsSlaves;
    MPI_Comm myComm;
    vector<int> vecSlave0;

    // === Initialize MPI and get ranke ===
    MPI_Init(&argc, &argv);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcessorsTotal);
    numProcessorsSlaves = numProcessorsTotal - 1;

    // === Correct the number of slaves if it is to heigh ===
    definitions::numSlaves = min(definitions::numSlaves, numProcessorsSlaves);

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
