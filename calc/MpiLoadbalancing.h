/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
    int size;
    int process;
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
        MPI_Comm mympicomm;
        MPI_Comm_split(mpi->world(),myColor,mpi->rank(), &mympicomm); // TODO integrate this into mympi
        mympi = std::make_shared<MpiVariables>(mympicomm);

        // === Construct new data type ===
        Triple trp;
        int cntBlock = 3; // 4
        int lenBlock[] = {1, 1, 1}; // 1, 1, 1, 1
        MPI_Aint displ[] = {0, reinterpret_cast<MPI_Aint>(&trp.col) - reinterpret_cast<MPI_Aint>(&trp), reinterpret_cast<MPI_Aint>(&trp.val) - reinterpret_cast<MPI_Aint>(&trp)}; // ..., sizeof(trp)
        MPI_Datatype type[] = {IDX_T, IDX_T, REAL_T}; //IDX_T, IDX_T, REAL_T, MPI_UB

        MPI_Type_create_struct(cntBlock, &lenBlock[0], &displ[0], &type[0], &MPI_COSTUME);
        MPI_Type_commit(&MPI_COSTUME);

        // === Wait for all datatypes and slave communicators to be build ===
        MPI_Barrier(mpi->world());
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

                // Do prework
                doPrework(numWork);

                // --- Increase number of working slaves0 ---
                numWorkingSlaves0 += 1;
            } else {
                // --- Receive message from slaves ---
                Message message_in = ReceiveFromSlave();
                size_t num = numSlave2numWork[message_in.process];

                // Save message
                dataOut[num]->deserialize(message_in.data);

                // Do postwork
                doPostwork(num);

                // --- Send new work ---
                numSlave2numWork[message_in.process]=numWork;
                Message message_out = Message(work, work.size(), message_in.process);
                SendToSlave(message_out);

                // Do prework
                doPrework(numWork);
            }
        }

        // === Collect - until now unreceived - results ===
        for (unsigned int i = 0; i < numWorkingSlaves0; i++) {
            // --- Receive message from slaves ---
            Message message_in = ReceiveFromSlave();
            size_t num = numSlave2numWork[message_in.process];

            // Save message
            dataOut[num]->deserialize(message_in.data);

            // Do rework
            doPostwork(num);
        }

        // === Kill all slaves0 ===
        std::vector<MPI_Request> requests(vecSlave0.size());

        size_t pos = 0;
        for (auto &dest: vecSlave0) {
            MPI_Isend(NULL, 0, MPI_INT, dest, DIETAG, mpi->world(), &requests[pos++]);
        }

        MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
    }

    virtual void runSlave() = 0;

    virtual void doPostwork(size_t numWork) {
        (void) numWork;
    }

    virtual void doPrework(size_t numWork) {
        (void) numWork;
    }

    void SendToSlave(Message &message) {
        MPI_Send(&message.data[0], message.size, BYTE_T, message.process, WORKTAG, mpi->world());
    }

    void SendToMaster(Message &message) {
        MPI_Send(&message.data[0], message.size, BYTE_T, 0, 0, mpi->world());
    }

    Message ReceiveFromSlave() {
        // Probe for an incoming message from slaves
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpi->world(), &status);

        // Initialize new message
        Message message;
        message.process = status.MPI_SOURCE;

        // Allocate buffer
        MPI_Get_count(&status,BYTE_T,&message.size);
        message.data.resize(message.size);

        // Receive message
        MPI_Recv(&message.data[0], message.size, BYTE_T, message.process, MPI_ANY_TAG, mpi->world(), &status);

        return message;
    }

    Message ReceiveFromMaster() {
        // Probe for an incoming message from master
        MPI_Status status;
        MPI_Probe(0, MPI_ANY_TAG, mpi->world(), &status);

        // Initialize new message
        Message message;
        message.process = 0;

        if (status.MPI_TAG != DIETAG) {
            // Allocate buffer
            MPI_Get_count(&status,BYTE_T,&message.size);
            message.data.resize(message.size);

            // Receive message
            MPI_Recv(&message.data[0], message.size, BYTE_T, 0, MPI_ANY_TAG, mpi->world(), &status);

        } else {
            // Receive message
            MPI_Recv(NULL, 0, BYTE_T, 0, MPI_ANY_TAG, mpi->world(), &status);

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
        MPI_Gather(&sz, 1, MPI_INT, &numElements[0], 1, MPI_INT, 0, mympi->world());

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
        MPI_Gatherv(&vectorIn[0], vectorIn.size(), MPI_COSTUME, &vectorOut[0], &numElements[0], &displacement[0], MPI_COSTUME, 0, mympi->world());

        return vectorOut;
    }

private:
    int numSlaves;
    std::vector<int> vecSlave0;
    MPI_Datatype MPI_COSTUME;
};

#endif
