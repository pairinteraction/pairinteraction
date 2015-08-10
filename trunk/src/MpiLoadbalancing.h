#ifndef MPILOADBALANCING_H
#define MPILOADBALANCING_H

#include "dtypes.h"

#include "MpiEnvironment.h"
#include "MpiVariables.h"

#include <math.h>
#include <vector>
#include <memory>

enum {WORKTAG = 0, DIETAG = 1};

struct Message{
    Message() : data(), size(0), process(-1){ }
    Message(bytes_t &data, size_t size, long process) : data(data), size(size), process(process){ }
    bytes_t data;
    long size;
    long process;
};

template <class TIn, class TOut>
class MpiLoadbalancing {
public:
    MpiLoadbalancing(std::shared_ptr<MpiEnvironment> mpi, int num) : mpi(mpi) {
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
        Triple trp;
        int cntBlock = 4;
        int lenBlock[] = {1, 1, 1, 1};
        MPI_Aint displ[] = {0, reinterpret_cast<MPI_Aint>(&trp.col) - reinterpret_cast<MPI_Aint>(&trp), reinterpret_cast<MPI_Aint>(&trp.val) - reinterpret_cast<MPI_Aint>(&trp), sizeof(trp)};
        MPI_Datatype type[] = {IDX_T, IDX_T, REAL_T, MPI_UB};

        MPI_Type_struct(cntBlock, &lenBlock[0], &displ[0], &type[0], &MPI_COSTUME);
        MPI_Type_commit(&MPI_COSTUME);

        // === Wait for all datatypes and slave communicators to be build ===
        mpi->world().Barrier();
    }

protected:
    std::shared_ptr<MpiVariables> mympi;
    std::shared_ptr<MpiEnvironment> mpi;

    void run(std::vector<std::shared_ptr<TIn> > &dataIn, std::vector<std::shared_ptr<TOut> > &dataOut) {
        if (this->mpi->rank() == 0) {
            this->runMaster(dataIn, dataOut);
        } else {
            this->runSlave();
        }
    }

    void runMaster(std::vector<std::shared_ptr<TIn> > &dataIn, std::vector<std::shared_ptr<TOut> > &dataOut) {
        size_t numWorkingSlaves0 = 0;
        std::vector<size_t> numSlave2numWork(mpi->size());

        dataOut.reserve(dataIn.size());
        for (size_t i = 0; i < dataIn.size(); ++i) {
            dataOut.push_back(std::move(std::make_shared<TOut>()));
        }

        // === Loop over work ===
        for (size_t numWork=0; numWork<dataIn.size(); ++numWork) {
            auto work = dataIn[numWork]->serialize();

            if (numWorkingSlaves0 < vecSlave0.size()) {
                // --- Send new work ---
                numSlave2numWork[vecSlave0[numWorkingSlaves0]]=numWork;
                Message message_out = Message(work, work.size(), vecSlave0[numWorkingSlaves0]);
                SendToSlave(message_out);

                // --- Increase number of working slaves0 ---
                numWorkingSlaves0 += 1;
            } else {
                // --- Receive message from slaves ---
                Message message_in = ReceiveFromSlave();

                // Save message
                dataOut[numSlave2numWork[message_in.process]]->deserialize(message_in.data);

                // --- Send new work ---
                numSlave2numWork[message_in.process]=numWork;
                Message message_out = Message(work, work.size(), message_in.process);
                SendToSlave(message_out);
            }
        }

        // === Collect - until now unreceived - results ===
        for (unsigned int i = 0; i < numWorkingSlaves0; i++) {
            // --- Receive message from slaves ---
            Message message_in = ReceiveFromSlave();

            // Save message
            dataOut[numSlave2numWork[message_in.process]]->deserialize(message_in.data);
        }

        // === Kill all slaves0 ===
        std::vector<MPI::Request> requests;
        requests.reserve(std::min(dataIn.size(), vecSlave0.size()));

        for (auto &dest: vecSlave0) {
            requests.push_back(mpi->world().Isend(NULL, 0, MPI::INT, dest, DIETAG));
        }

        MPI::Request::Waitall(requests.size(), &requests[0]);
    }

    virtual void runSlave() = 0;

    void SendToSlave(Message &message) {
        mpi->world().Send(&message.data[0], message.size, BYTE_T, message.process, WORKTAG);
    }

    void SendToMaster(Message &message) {
        mpi->world().Send(&message.data[0], message.size, BYTE_T, 0, 0);
    }

    Message ReceiveFromSlave() {
        // Probe for an incoming message from slaves
        MPI::Status status;
        mpi->world().Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, status);

        // Initialize new message
        Message message;
        message.process = status.Get_source();

        // Allocate buffer
        message.size = status.Get_count(BYTE_T);
        message.data.resize(message.size);

        // Receive message
        mpi->world().Recv(&message.data[0], message.size, BYTE_T, message.process, MPI::ANY_TAG);

        return message;
    }

    Message ReceiveFromMaster() {
        // Probe for an incoming message from master
        MPI::Status status;
        mpi->world().Probe(0, MPI::ANY_TAG, status);

        // Initialize new message
        Message message;
        message.process = 0;

        if (status.Get_tag() != DIETAG) {
            // Allocate buffer
            message.size = status.Get_count(BYTE_T);
            message.data.resize(message.size);

            // Receive message
            this->mpi->world().Recv(&message.data[0], message.size, BYTE_T, 0, MPI::ANY_TAG);

        } else {
            // Receive message
            this->mpi->world().Recv(NULL, 0, BYTE_T, 0, MPI::ANY_TAG);

            // Account for the DIETAG
            message.size = -1;
        }

        return message;
    }

    std::vector<Triple> Gathervector(std::vector<Triple> &vectorIn) {
        std::vector<Triple> vectorOut;

        // Tell slave0 how many elements one has
        std::vector<int> numElements(mympi->size());
        int sz = vectorIn.size();
        mympi->world().Gather(&sz, 1, MPI::INT, &numElements[0], 1, MPI::INT, 0);

        // Calculate displacements
        std::vector<int> displacement(mympi->size());
        if ( mympi->rank() == 0 ) {

            int nTotal = 0;
            for (int i = 0; i < mympi->size(); i++) {
                displacement[i] = nTotal;
                nTotal += numElements[i];
            }

            // Allocate buffer
            vectorOut.resize(nTotal);
        }

        // Collect everything into slave0
        mympi->world().Gatherv(&vectorIn[0], vectorIn.size(), MPI_COSTUME, &vectorOut[0], &numElements[0], &displacement[0], MPI_COSTUME, 0);

        return vectorOut;
    }

private:
    int numSlaves;
    std::vector<int> vecSlave0;
    MPI_Datatype MPI_COSTUME;
};

#endif
