#include "MpiLoadbalancing.h"

MpiLoadbalancing::MpiLoadbalancing(std::shared_ptr<MpiEnvironment> mpi, int num) : mpi(mpi) {
    // === Correct the number of slaves if it is to heigh ===
    numSlaves = std::min(num, mpi->size()-1);

    // === Construct a vector that contains the rank of the 0th slave inside each slave group ===
    int numUnusedProcessors = mpi->size()-1;
    int numUnusedSlaves = numSlaves;
    int idxProcessor = 1;

    int myColor = numSlaves;
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
    mympi = std::make_shared<MpiVariables>(mpi->world().Split(myColor, mpi->rank()));

    // === Construct new data type ===
    triple trp;
    int cntBlock = 4;
    int lenBlock[cntBlock] = {1, 1, 1, 1};
    MPI_Aint displ[cntBlock] = {0, reinterpret_cast<MPI_Aint>(&trp.col) - reinterpret_cast<MPI_Aint>(&trp), reinterpret_cast<MPI_Aint>(&trp.val) - reinterpret_cast<MPI_Aint>(&trp), sizeof(trp)};
    MPI_Datatype type[cntBlock] = {DIDX, DIDX, DREAL, MPI_UB};

    MPI_Type_struct(cntBlock, &lenBlock[0], &displ[0], &type[0], &MPI_COSTUME);
    MPI_Type_commit(&MPI_COSTUME);

    // === Wait for all datatypes and slave communicators to be build ===
    mpi->world().Barrier();
}

void MpiLoadbalancing::runMaster(std::vector<std::shared_ptr<Serializable> > &dataIn, std::vector<std::shared_ptr<Serializable> > &dataOut) {
    size_t numWorkingSlaves0 = 0;
    std::vector<size_t> numSlave2numWork(mpi->size());

    // === Loop over work ===
    for (size_t numWork=0; numWork<dataIn.size(); ++numWork) {
        auto work = dataIn[numWork]->serialize();

        if (numWorkingSlaves0 < vecSlave0.size()) {
            // --- Send new work ---
            int numSlave = vecSlave0[numWorkingSlaves0];
            numSlave2numWork[numSlave]=numWork;

            mpi->world().Send(&work[0], work.size(), MPI::UNSIGNED_CHAR, numSlave, WORKTAG);

            // --- Increase number of working slaves0 ---
            numWorkingSlaves0 += 1;
        } else {
            // --- Receive message from slaves ---

            // Probe for an incoming message from master
            MPI::Status status;
            mpi->world().Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, status);

            // Allocate buffer
            int msgSize = status.Get_count(MPI::UNSIGNED_CHAR);
            std::vector<unsigned char> messagebuffer(msgSize);

            // Receive message
            mpi->world().Recv(&messagebuffer[0], msgSize, MPI::UNSIGNED_CHAR, status.Get_source(), MPI::ANY_TAG);

            // Save message
            dataOut[numSlave2numWork[status.Get_source()]]->deserialize(messagebuffer);

            // --- Send new work ---
            int numSlave = status.Get_source();
            numSlave2numWork[numSlave]=numWork;

            mpi->world().Send(&work[0], work.size(), MPI::UNSIGNED_CHAR, numSlave, WORKTAG);
        }
    }

    // === Collect - until now unreceived - results ===
    for (unsigned int i = 0; i < numWorkingSlaves0; i++) {
        // --- Receive message from slaves ---

        // Probe for an incoming message from master
        MPI::Status status;
        mpi->world().Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, status);

        // Allocate buffer
        int msgSize = status.Get_count(MPI::UNSIGNED_CHAR);
        std::vector<unsigned char> messagebuffer(msgSize);

        // Receive message
        mpi->world().Recv(&messagebuffer[0], msgSize, MPI::UNSIGNED_CHAR, status.Get_source(), MPI::ANY_TAG);

        // Save message
        dataOut[numSlave2numWork[status.Get_source()]]->deserialize(messagebuffer);
    }

    // === Kill all slaves0 ===
    std::vector<MPI::Request> requests;
    requests.reserve(std::min(dataIn.size(), vecSlave0.size()));

    for (auto &dest: vecSlave0) {
        requests.push_back(mpi->world().Isend(NULL, 0, MPI::INT, dest, DIETAG));
    }

    MPI::Request::Waitall(requests.size(), &requests[0]);
}

void MpiLoadbalancing::runSlave(std::shared_ptr<Serializable> bufferSerializable, std::shared_ptr<Vectorizable> bufferVectorizable) {
    // === Do work ===
    while (true) {
        std::vector<unsigned char> work;
        int msgSize;

        // --- Receive message from master ---
        if (mympi->rank() == 0) {
            // Probe for an incoming message from master
            MPI::Status status;
            mpi->world().Probe(0, MPI::ANY_TAG, status);

            if (status.Get_tag() != DIETAG) {
                // Get message size
                msgSize = status.Get_count(MPI::UNSIGNED_CHAR);

                // Allocate buffer
                work.resize(msgSize);

                // Receive message
                mpi->world().Recv(&work[0], work.size(), MPI::UNSIGNED_CHAR, 0, MPI::ANY_TAG);
            } else {
                // Receive message
                mpi->world().Recv(NULL, 0, MPI::UNSIGNED_CHAR, 0, MPI::ANY_TAG);

                // Set message size to -1 to account for the DIETAG
                msgSize = -1;
            }
        }

        // --- Receive message from slave0 ---
        if (mympi->size() > 1) {
            // Broadcast message size
            mympi->world().Bcast(&msgSize, 1, MPI::INT, 0);

            if (msgSize != -1) {
                if ( mympi->rank() != 0 ) {
                    // Allocate buffer
                    work.resize(msgSize);
                }

                // Broadcast message
                mympi->world().Bcast(&work[0], work.size(), MPI::UNSIGNED_CHAR, 0);
            }
        }


        // --- Stop working in case of DIETAG ---
        if (msgSize == -1) {

            break;
        }

        // --- Do calculations ---
        std::cout << "Hello from rank " << mpi->rank() <<  ", rank within slave group " << mympi->rank() << std::endl;
        bufferSerializable->deserialize(work);
        auto resultunvectorized = doMainprocessing(bufferSerializable);

        // --- Send results to slave0 ---
        std::vector<triple> collection;

        if (resultunvectorized != NULL) {
            auto result = resultunvectorized->vectorize();

            if (mympi->size() > 1) {
                // Tell slave0 how many elements one has
                std::vector<int> numElements(mympi->size());
                int sz = result.size();
                mympi->world().Gather(&sz, 1, MPI::INT, &numElements[0], 1, MPI::INT, 0);

                std::vector<int> displacement(mympi->size());
                if ( mympi->rank() == 0 ) {
                    // Calculate displacements
                    int nTotal = 0;
                    for (int i = 0; i < mympi->size(); i++) {
                        displacement[i] = nTotal;
                        nTotal += numElements[i];
                    }

                    // Allocate buffer
                    collection.resize(nTotal);
                }

                // Collect everything into slave0
                mympi->world().Gatherv(&result[0], result.size(), MPI_COSTUME, &collection[0], &numElements[0], &displacement[0], MPI_COSTUME, 0);
            } else {
                collection = result;
            }
        } else {
            bufferVectorizable = NULL;
        }

        // --- Do postprocessing and send results to master ---
        if (mympi->rank() == 0) {
            // Do postprocessing
            if (bufferVectorizable != NULL) {
                bufferVectorizable->devectorize(collection);
            }
            auto results = doPostprocessing(bufferVectorizable, bufferSerializable)->serialize();

            // Send results to master
            mpi->world().Send(&results[0], results.size(),  MPI::UNSIGNED_CHAR, 0, 0);
        }
    }
}


void MpiLoadbalancing::run(std::vector<std::shared_ptr<Serializable>> &dataIn, std::vector<std::shared_ptr<Serializable>> &dataOut, std::shared_ptr<Serializable> bufferSerializable, std::shared_ptr<Vectorizable> bufferVectorizable) {
    if (mpi->rank() == 0) {
        this->runMaster(dataIn, dataOut);
    } else {
        this->runSlave(bufferSerializable, bufferVectorizable);
    }
}
